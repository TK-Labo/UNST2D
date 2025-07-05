# UNST2D

## Introduction
このプロジェクトでは、非構造格子二次元不定流モデル (UNST2D: Unstructured grid 2D unsteady flow model) に、下水道・圃場整備エリア、田んぼダム、樹林帯（防備林）、連続盛土、一次元河道などの要素モデルを追加し、流域治水対策の効果を見える化するためのモデルを開発しています。国立研究開発法人土木研究所が公開する降雨流出氾濫（RRI）モデルのほか分布型流出モデル等の出力を計算領域の外縁に与えることも可能です。

## Citation
このコードを利用した計算結果の公表・頒布に際しては、以下の論文を引用してください。  
Please cite the following paper when publishing or distributing calculation results using this code.

川池 健司, 井上 和也, 戸田 圭一 (2000) 非構造格子の都市氾濫解析への適用, 水工学論文集, 44: 461-466.  
Kenji KAWAIKE, Kazuya INOUE, Kei-ichi TODA (2000) Applications of unstructured meshes to inundation flow analysis in urban area, Annual journal of Hydraulic Engineering, JSCE, 44, 461-466.  
https://doi.org/10.2208/prohe.44.461

山村 孝輝, 瀧健太郎 (2024) 遺伝的アルゴリズムを用いた田んぼダムの最適配置探索法の提案，土木学会論文集, 81(16): 24-16046.  
Koki YAMAMURA, Shunji NISHINO, Masafumi YAMADA, Takahiro SAYAMA, Kenji KAWAIKE, Kentaro TAKI (2025) PROPOSAL OF SEARCH METHOD FOR OPTIMAL PLACEMENT OF “PADDY-FIELD-DAM” USING GENETIC ALGORITHM, Japanese Journal of JSCE, 31, 81(16): 24-16046.  
https://doi.org/10.2208/jscejj.24-16046  
  
## Contents

プロジェクトは次の主要なディレクトリとファイルで構成されています:

- `src/`: Fortranソースコードが格納されたディレクトリ
  - UNSTモデル関連のファイル（`UNST*.f90`）
  - Makefile

## Compile
UNST2Dによる計算を行うプログラム`UNST.exe`を生成:  

```bash
make
```
  
## Run
必要な入力ファイルを準備した後、以下のコマンドを実行すると計算が始まります:
```bash
./RRI_UNST.exe
```
  
## License

Copyright (c) 2025 K.Kawaike & TK Labo

Licensed under the [MIT](https://github.com/TK-Labo/RRI2UNST2D/blob/main/LICENSE) license.

## Coded by

Kenji Kawaike & TK Labo

TK Labo Members:

- Koki Yamamura 2021-2025
- Daiki Baba 2022-
- Shunji Nishino 2023-
- Kentaro Taki 2017-
