# UNST2D

## Introduction
このプロジェクトでは、非構造格子二次元不定流モデル (UNST2D: Unstructured grid 2D unsteady flow model) に、下水道・圃場整備エリア、田んぼダム、樹林帯（防備林）、連続盛土、一次元河道などの要素モデルを追加し、流域治水対策の効果を見える化するためのモデルを開発しています。国立研究開発法人土木研究所が公開する降雨流出氾濫（RRI）モデルのほか分布型流出モデル等の出力を計算領域の外縁に与えることも可能です。

In this project, we are developing a model to visualize the effects of river basin flood control measures by adding element models such as sewerage and farmland improvement areas, rice paddy dams, forest belts (defensive forests), continuous embankments, and one-dimensional river channels to the Unstructured Grid 2D Unsteady Flow Model (UNST2D). It is also possible to provide the output of the Rainfall-Runoff-Inundation (RRI) model released by the Public Works Research Institute, as well as distributed runoff models, to the outer edge of the calculation domain.

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
  - UNST2Dモデル関連のファイル（`UNST*.f90`）
  - Makefile

The project consists of the following main directories and files:

- `src/`: Directory where Fortran source code is stored
  - UNST2D model related files（`UNST*.f90`）
  - Makefile

## Compile
UNST2Dによる計算を行うプログラム`UNST.exe`を生成:  
Generate the program `UNST.exe` to perform calculations using UNST2D:  

```bash
make
```
  
## Run  
必要な入力ファイルを準備した後、以下のコマンドを実行すると計算が始まります:
After preparing the necessary input files, run the following command to start the calculation:
```bash
./RRI_UNST.exe
```
  
## License

Licensed under the [MIT](https://github.com/TK-Labo/RRI2UNST2D/blob/main/LICENSE) license.

Copyright (c) 2025 K.Kawaike & TK Labo

## Coded by

Kenji Kawaike & TK Labo

TK Labo Members:

- Koki Yamamura 2021-2025
- Daiki Baba 2022-
- Shunji Nishino 2023-
- Kentaro Taki 2017-
