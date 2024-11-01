import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os

class Scanpy():
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "file_path": ("STRING", {"multiline": False, "default": ""}),
            }
        }
    
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("result", "ANNDATA")
    FUNCTION = "Data_on"
    CATEGORY = "Bioinformatics Nodes"
    
    def Data_on(self, file_path):
        sc.settings.set_figure_params(dpi=50, facecolor="white")

        adatas = {}

        # file_path를 사용하여 파일 경로 사전 생성
        file_paths = {
            "s1d1": f"{file_path}/s1d1_filtered_feature_bc_matrix.h5",
            "s1d3": f"{file_path}/s1d3_filtered_feature_bc_matrix.h5",
        }

        for sample_id, path in file_paths.items():
            sample_adata = sc.read_10x_h5(path)
            sample_adata.var_names_make_unique()
            adatas[sample_id] = sample_adata
            
        # 모든 샘플 데이터를 하나의 AnnData 객체로 결합
        adata = ad.concat(adatas, label="sample")
        adata.obs_names_make_unique()

        # 미토콘드리아, 리보솜, 헤모글로빈 유전자 그룹을 정의
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

        # 품질 관리 메트릭스 계산
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
        )
        
        # 결과 반환
        result_counts = adata.obs["sample"].value_counts().to_string()
        return (result_counts, adata)
    
class Plot_Violin():
    def __init__(self):
        pass
    
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
                "save_path": ("STRING", {"multiline": False, "default": ""}),
                "file_name": ("STRING", {"multiline": False, "default": "violin"})
            },
        }
    
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("Path",)
    FUNCTION = "plot_violin"
    CATEGORY = "Bioinformatics Nodes"
    
    def plot_violin(self, adata, save_path, file_name):
        # 지정된 경로가 존재하지 않으면 생성
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        
        # 바이올린 플롯 생성
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False  # show를 False로 설정하여 바로 표시하지 않음
        )

        # save_path와 file_name을 결합하여 파일을 저장
        full_path = os.path.join(save_path, f"{file_name}.png")  # 확장자 추가
        plt.savefig(full_path)
        plt.close()  # 플롯을 닫아 메모리 절약

        return (f"{full_path}",)
    
class Plot_Scatter():
    def __init__(self):
        pass
    
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
                "save_path": ("STRING", {"multiline": False, "default": ""}),
                "file_name": ("STRING", {"multiline": False, "default": "scatter"})
            },
        }
    
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("Path",)
    FUNCTION = "plot_scatter"
    CATEGORY = "Bioinformatics Nodes"
    
    def plot_scatter(self, adata, save_path, file_name):
        # 지정된 경로가 존재하지 않으면 생성
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        
        # 산점도 생성
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="n_genes_by_counts",
            color="pct_counts_mt",
            show=False  # 바로 표시하지 않음
        )

        # save_path와 file_name을 결합하여 파일을 저장
        full_path = os.path.join(save_path, f"{file_name}.png")  # 확장자 추가
        plt.savefig(full_path)
        plt.close()  # 플롯을 닫아 메모리 절약

        return (f"{full_path}",)

class FilterData:
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
                "min_genes": ("INT", {"default": 100}),
                "min_cells": ("INT", {"default": 3}),
            },
        }

    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("ANNDATA",)
    FUNCTION = "filter_data"
    CATEGORY = "Bioinformatics Nodes"
    
    def filter_data(self, adata, min_genes, min_cells):
        # 세포 필터링: 최소 min_genes 수의 유전자가 발현된 세포만 남김
        sc.pp.filter_cells(adata, min_genes=min_genes)
        
        # 유전자 필터링: 최소 min_cells 수의 세포에서 발현된 유전자만 남김
        sc.pp.filter_genes(adata, min_cells=min_cells)
        
        return (adata,)

class DetectDoublets:
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
                "batch_key": ("STRING", {"default": "sample"})  # 기본 값은 "sample"
            }
        }

    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("result", "ANNDATA")
    FUNCTION = "detect_doublets"
    CATEGORY = "Bioinformatics Nodes"

    def detect_doublets(self, adata, batch_key):
        # Scrublet을 사용하여 이중 세포 감지
        sc.pp.scrublet(adata, batch_key=batch_key)

        # 결과 반환 메시지
        return adata, f"Doublet detection completed using batch_key='{batch_key}'"

class PreprocessRNASeq:
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
            },
        }

    RETURN_TYPES = ("STRING","STRING")
    RETURN_NAMES = ("result", "ANNDATA")
    FUNCTION = "preprocess_data"
    CATEGORY = "Bioinformatics Nodes"

    def preprocess_data(self, adata):
        # Saving count data
        adata.layers["counts"] = adata.X.copy()

        # Normalizing to median total counts
        sc.pp.normalize_total(adata)

        # Logarithmize the data
        sc.pp.log1p(adata)

        return (adata, "Preprocessing completed.",)

class IdentifyAndPlotHighlyVariableGenes:
    def __init__(self):
        pass

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
                "n_top_genes": ("INT", {"default": 2000}),
                "batch_key": ("STRING", {"default": "sample"}),
                "save_path": ("STRING", {"multiline": False, "default": ""}),
                "file_name": ("STRING", {"multiline": False, "default": "highly_variable_genes"})
            },
        }

    RETURN_TYPES = ("STRING","STRING")
    RETURN_NAMES = ("Path", "ANNDATA")
    FUNCTION = "identify_and_plot_hvg"
    CATEGORY = "Bioinformatics Nodes"

    def identify_and_plot_hvg(self, adata, n_top_genes, batch_key, save_path, file_name):
        # 고변동 유전자 선택
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)

        # 지정된 경로가 존재하지 않으면 생성
        folder_path = os.path.dirname(save_path)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        # 고변동 유전자 시각화
        sc.pl.highly_variable_genes(adata, show=False)

        # matplotlib의 savefig를 사용하여 지정된 경로에 직접 저장
        plt.savefig(os.path.join(save_path, f"{file_name}.png"))
        plt.close()  # 플롯을 닫아 메모리 절약

        return (adata, f"{os.path.join(save_path, f'{file_name}.png')}",)
#adata_hvg / adata 나누기
class PCAAnalysis:
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
                "n_top_genes": ("INT", {"default": 2000}),
                "n_pcs": ("INT", {"default": 50}),
                "save_path": ("STRING", {"multiline": False, "default": "output"}),
                "file_name": ("STRING", {"multiline": False, "default": "scatter"}),
            },
        }
    
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("ANNDATA_HVG", "ANNDATA",)
    FUNCTION = "perform_pca_analysis"
    CATEGORY = "Bioinformatics Nodes"
    
    def perform_pca_analysis(self, adata, n_top_genes, n_pcs, save_path, file_name):
        # 고변동 유전자 식별
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key="sample")
        
        # 고변동 유전자만으로 adata 업데이트 (실제 데이터로 변환)
        adata_hvg = adata[:, adata.var['highly_variable']].copy()
        
        # PCA 수행
        sc.tl.pca(adata_hvg, n_comps=n_pcs)
        
        # PCA 분산 비율 플롯 생성
        sc.pl.pca_variance_ratio(adata_hvg, n_pcs=n_pcs, log=True, show=False)

        # 지정된 경로에 저장할 폴더가 없으면 생성
        folder_path = os.path.dirname(save_path)
        if folder_path and not os.path.exists(folder_path):
            os.makedirs(folder_path)

        full_path = os.path.join(save_path, f"{file_name}.png")  # 확장자 추가
        plt.savefig(full_path)
        plt.close()
        
        return (adata_hvg, adata)

class PCAPlot:
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "adata": ("STRING", {"forceInput": True}),
                "adata_hvg": ("STRING", {"forceInput": True}),
                "color": ("STRING", {"multiline": True, "default": "sample, sample, pct_counts_mt, pct_counts_mt"}),
                "dimensions": ("STRING", {"multiline": True, "default": "0,1;2,3"}),
                "ncols": ("INT", {"default": 2}),
                "size": ("FLOAT", {"default": 2.0}),
                "save_path": ("STRING", {"multiline": False, "default": "output"}),
                "file_name": ("STRING", {"multiline": False, "default": "file name"}),
            },
        }

    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("ANNDATA",)
    FUNCTION = "plot_pca"
    CATEGORY = "Bioinformatics Nodes"
    
    def plot_pca(self, adata, adata_hvg, color, dimensions, ncols, size, save_path, file_name):
        color_list = color.split(",")  # 문자열을 리스트로 변환
        dimension_pairs = [tuple(map(int, dim.split(","))) for dim in dimensions.split(";")]  # 문자열을 차원 쌍으로 변환
        if adata_hvg : 
            # PCA 플롯 생성
            sc.pl.pca(
                adata_hvg,
                color=color_list,
                dimensions=dimension_pairs,
                ncols=ncols,
                size=size,
                show=False  # 바로 보여주지 않도록 설정
            )

            # 지정된 경로에 저장할 폴더가 없으면 생성
            folder_path = os.path.dirname(save_path)
            if folder_path and not os.path.exists(folder_path):
                os.makedirs(folder_path)

            full_path = os.path.join(save_path, f"{file_name}.png")
            plt.savefig(full_path)
            plt.close()

            return f"PCA plot saved at: {save_path}"
        
        else:
            # PCA 플롯 생성
            sc.pl.pca(
                adata,
                color=color_list,
                dimensions=dimension_pairs,
                ncols=ncols,
                size=size,
                show=False  # 바로 보여주지 않도록 설정
            )

            # 지정된 경로에 저장할 폴더가 없으면 생성
            folder_path = os.path.dirname(save_path)
            if folder_path and not os.path.exists(folder_path):
                os.makedirs(folder_path)

            full_path = os.path.join(save_path, f"{file_name}.png")
            plt.savefig(full_path)
            plt.close()

            return f"PCA plot saved at: {save_path}"
    
NODE_CLASS_MAPPINGS = {
    "Scanpy" : Scanpy,
    "Plot_Violin" : Plot_Violin,
    "Plot_Scatter" : Plot_Scatter,
    "FilterData" : FilterData,
    "DetectDoublets": DetectDoublets,
    "PreprocessRNASeq": PreprocessRNASeq,
    "IdentifyAndPlotHighlyVariableGenes": IdentifyAndPlotHighlyVariableGenes,
    "PCAAnalysis": PCAAnalysis,
    "PCAPlot": PCAPlot,
    }

NODE_DISPLAY_NAME_MAPPINGS = {
    "Scanpy" : "Scanpy",
    "Plot_Violin" : "Plot Violin",
    "Plot_Scatter" : "Plot Scatter",
    "FilterData" : "Filter Data",
    "DetectDoublets" : "Detect Doublets",
    "PreprocessRNASeq" : "Preprocess RNASeq",
    "IdentifyAndPlotHighlyVariableGenes" : "Plot Identify Highly Variable Genes",
    "PCAAnalysis" : "PCA Analysis",
    "PCAPlot": "PCAPlot",
    }

#plot관련은 통일시켜서 하나의 노드에서 선택에 따라 다른 plot그리게끔하기