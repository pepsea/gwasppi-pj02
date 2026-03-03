# GWAS-PPI Analysis Pipeline (pj02)

## 概要
このプロジェクト（`pj02`）は、**「疾患はゲノム構造の影響を受けた下流の遺伝子が引き起こす生体メカニズムの異常」** という仮説に基づき、GWAS (ゲノムワイド関連解析) と PPI (タンパク質間相互作用) ネットワークを統合し、疾患の上流原因（ゲノム）から下流結果（メカニズム異常）の中核を担う重要な遺伝子群を網羅的に抽出し、新規の創薬ターゲットを探索するための統合解析パイプラインです。

## 全体設計とデータフロー

本システムは以下のステップで構成されています。これらは Jupyter Notebook (`notebooks/gwas_ppi_analysis.ipynb`) から一元的に実行可能です。

### 1. GWAS SNP → 遺伝子マッピング (`gwas_fetcher.py`)
- **目的:** 特定疾患・形質のリスクとして報告された SNP（一塩基多型）を取得し、関連する遺伝子群を特定する。
- **データソース:** [GWAS Catalog REST API](https://www.ebi.ac.uk/gwas/rest/docs/api)
- **処理内容:**
  - 疾患名（例: "Type 2 diabetes"）または EFO ID でクエリを実行。
  - 大量のアソシエーションデータから上位のスタディ（`max_studies`パラメータで制限可能）に含まれる有意な SNP（通常 P値 < 5e-8）を抽出。
  - 各 SNP の P値（$-\log_{10}P$）を保持し、周辺の Mapped Gene を取得・マッピング。

### 2. PPI ネットワークの構築 (`ppi_fetcher.py`)
- **目的:** GWASで特定した遺伝子群（Seed Genes）から、下流に繋がる生体メカニズムのネットワークを構築する。
- **データソース:** OmniPathDB (SIGNOR, Reactome), STRING API
- **処理内容:**
  - Seed Genes を起点に、相互作用が報告されている遺伝子をN階層（デフォルト1～3）まで探索。
  - 複数のソースからPPIを取得し、エッジごとに「ソース」をメタデータとして保持（STRINGの場合はスコアによる足切りも可能）。

### 3. バリアント機能アノテーション (`gene_scorer.py`)
- **目的:** SNPが単なるマーカーではなく、タンパク質機能に対して直接的な変化（Loss of Function / Gain of Function）をもたらすかを評価する。
- **データソース:** [Ensembl VEP API](https://rest.ensembl.org/documentation/info/vep_id_post)
- **処理内容:**
  - 各 SNP に対し VEP (Variant Effect Predictor) を叩き、転写産物に対する影響（Consequence）を取得。
  - `stop_gained`, `frameshift_variant` は **LoF**。
  - `missense_variant` で特定のスコア（SIFT, PolyPhen）が有害・病的と予測されるものは **GoF/Missense**。
  - それ以外を Coding / Regulatory などに分類。

### 4. フローベースの重要遺伝子スコアリング (`network_analysis.py`)
- **目的:** ネットワークの中心にある単なる「ハブ遺伝子」だけでなく、GWASからの「太い流れ（フロー）」の中にいる遺伝子を高く評価する。
- **処理内容:**
  - **初期流量の定義:** 各 GWAS 遺伝子に対して、 GWAS P値の対数スケールと、VEPで得たバリアント機能スコア（LoF等には高い重み）の合計値を初期スコア（Source Flow）として設定。
  - **フローの伝播 (Flow Propagation):** ランダムウォークの応用（PageRank的アプローチ）により、初期流量を持つ GWAS遺伝子から PPI エッジを通じて隣接する遺伝子へとスコアを伝播（減衰係数を伴い数ステップ拡散）。
  - **統合スコアリング:** 最終的な `flow_score` (ウェイト比重 大) に加え、ネットワークの次数・媒介中心性などの各種中心性指標を組み合わせ（Integrated Score）、スコア上位の遺伝子を「Key Genes」として抽出。

### 5. エンリッチメント解析 (`enrichment.py`)
- **目的:** 抽出された Key Genes が集積している生物学的パスウェイや疾患メカニズムを明らかにする。
- **データソース:** ローカル GMT ファイル (Reactome, GO, HPO)
- **処理内容:**
  - `gseapy` を用いずに、ローカルに保持した GMT ファイル群に対して Fisher's Exact Test を直接実行。
  - Benjamini-Hochberg (FDR) により多重検定補正を行い、有意なタームを抽出。

### 6. 創薬ターゲットの取得・オーバーレイ (`target_fetcher.py`)
- **目的:** ネットワークとパスウェイの中に、すでに開発されている、または未開拓のターゲット遺伝子がないかを照合する。
- **データソース:** [Pharos GraphQL API](https://pharos.nih.gov/), [ChEMBL API](https://www.ebi.ac.uk/chembl/api/data/docs)
- **処理内容:**
  - EFO ID または MeSH ID を用いて、当該疾患の適応症を持つ承認薬（Ph 4）および臨床試験中（Ph 1-3）の化合物と、その作用機序（ターゲット遺伝子）を API で動的に取得。
  - 双方のデータベースの情報をマージし、最大開発フェーズ（Max Phase）情報を生成。取得したリストは可視化モジュールに渡される。

### 7. 可視化 (`visualization.py`)
- **ツール:** Plotly
- **出力物:** 全て `.html` （インタラクティブビュー）で `output/` に出力。
  - **Sankey Diagram:** GWAS SNP → GWAS遺伝子 → Key PPI遺伝子 → パスウェイ への「情報の流れ」を可視化。
  - **PPI Network Plot:** スプリングレイアウトを用いたネットワーク図。ノードの大きさと色はスコアに依存。**創薬ターゲット遺伝子（💊）はマゼンタカラーの太線でハイライト**される。
  - **Pathway Network Plots:** 有意なパスウェイごとに切り出した部分ネットワーク群。
  - **Pathway-Gene Matrix:** 上位パスウェイ × 重要遺伝子群の参加状況を示すヒートマップ。ターゲットは `💊GENE (Ph4)` 等の形式で注釈される。
  - **Bubble / Bar Charts:** エンリッチメント結果の可視化。

---

## ディレクトリ構成
```text
pj02/
├── config.py                # グローバル設定・API URL・定数パラメータ (max_studies等)
├── id_mapper.py             # 遺伝子シンボル/Entrez ID間の相互変換
├── gwas_fetcher.py          # GWAS CatalogからのSNP取得機能
├── ppi_fetcher.py           # OmniPath/STRINGからの相互作用取得機能
├── gene_scorer.py           # Ensembl VEP連携・バリアントベースのスコアリング
├── target_fetcher.py        # Pharos & ChEMBL: EFO/MeSHからのターゲット取得
├── network_analysis.py      # フロー伝播・ネットワーク中心性解析
├── enrichment.py            # Local Fisher's Exact Test 解析
├── visualization.py         # 各種インタラクティブ可視化ビューの生成
├── notebooks/
│   └── gwas_ppi_analysis.ipynb # メイン実行パイプライン
├── data/                    # IDマッピングやGMT等の内部データ
└── output/                  # 解析の生成物 (CSV, HTML) が保存される先
```

## 環境構築・実行方法
詳細は既存のドキュメントに譲りますが、本一連の解析は `notebooks/gwas_ppi_analysis.ipynb` を上から下まで実行することで完結します。

各種パラメータ（探索層の深さ、対象疾患 `EFO_TRAIT_ID`、GWASスタディの取得制限数 `DEFAULT_MAX_GWAS_STUDIES`、スコアリングのウェイト等）は、Notebookの先頭セル、または `config.py` 内で一元的に調整可能です。
