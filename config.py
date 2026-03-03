"""
GWAS-PPI Pipeline Configuration (pj02)
=======================================
API キー、URL、デフォルトパラメータを管理
"""

import os

# ============================================================
# ディレクトリ設定
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# GWAS Catalog API
# ============================================================
GWAS_CATALOG_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"
GWAS_P_VALUE_THRESHOLD = 5e-8  # genome-wide significance

# ============================================================
# PPI Databases
# ============================================================

# SIGNOR
SIGNOR_DOWNLOAD_URL = "https://signor.uniroma2.it/download_entity.php"
SIGNOR_ALL_DATA_URL = "https://signor.uniroma2.it/Scripts/download_entity_csv.php"

# BioGRID
BIOGRID_BASE_URL = "https://webservice.thebiogrid.org/interactions/"
BIOGRID_ACCESS_KEY = os.environ.get("BIOGRID_ACCESS_KEY", "e8f7cf92fae447b3e9d6729aa6a815eb")

# STRING
STRING_BASE_URL = "https://string-db.org/api"
STRING_SPECIES = 9606
STRING_MIN_SCORE = 700  # デフォルト: Score > 0.7 (高信頼)

# Reactome
REACTOME_CONTENT_URL = "https://reactome.org/ContentService"
REACTOME_ANALYSIS_URL = "https://reactome.org/AnalysisService"

# ============================================================
# ID Mapping
# ============================================================
UNIPROT_ID_MAPPING_URL = "https://rest.uniprot.org/idmapping"

# ============================================================
# Ensembl VEP
# ============================================================
ENSEMBL_REST_URL = "https://rest.ensembl.org"

# ============================================================
# Enrichment GMTデータ
# ============================================================
REACTOME_GMT_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
GO_GMT_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c5.go.v2024.1.Hs.symbols.gmt"
HPO_ANNOTATION_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt"

# ============================================================
# デフォルト設定
# ============================================================
DEFAULT_ORGANISM = "Homo sapiens"
DEFAULT_TAXID = 9606
ENRICHMENT_FDR_THRESHOLD = 0.05
BACKGROUND_GENE_COUNT = 20000

# GWAS Studyの最大取得数 (P値が小さい = 有意なものから上位N件に絞る)
DEFAULT_MAX_GWAS_STUDIES = 50

# PPI 階層 (デフォルト: 1階層)
DEFAULT_PPI_LAYERS = 1

# デフォルト PPI ソース
DEFAULT_PPI_SOURCES = ["signor", "reactome", "string"]

# スコアリング重み
SCORING_WEIGHTS = {
    "gwas_pvalue": 1.0,       # -log10(p_value) を基本スコア (上限なし)
    "lof": 40.0,              # Loss of Function (P値 1e-40 相当の影響力)
    "gof": 35.0,              # Gain of Function
    "missense": 20.0,         # ミスセンスバリアント
    "coding_variant": 15.0,   # コーディング領域バリアント
    "regulatory": 10.0,       # 調節領域バリアント
    "ppi_multi_source": 10.0, # PPI 複数ソース確認 (1ソースあたり +10)
}
