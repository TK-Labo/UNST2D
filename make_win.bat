@echo off
setlocal enabledelayedexpansion
REM==============================================================================
REM RRI-UNST連携モデル用make.bat (Windows)
REM
REM 使い方:
REM   .\make_win.bat        : プログラムをコンパイルして実行ファイルを生成
REM
REM 最終更新:
REM==============================================================================
REM -------------------------------------------
REM  コンパイラの選択
REM -------------------------------------------
REM set COMPILER=ifx
set COMPILER=gfortran

REM コンパイラごとの設定
if /I "%COMPILER%"=="gfortran" (
    set FC=gfortran
    set OPTIONS=-O2 -fopenmp -cpp -ffree-line-length-none
    set MOD_FLAG=-J
    set OBJ_FLAG=-o
    set OSUFFIX=.o
) else (
    REM Intel Fortran (ifx)
    set FC=ifx
    set OPTIONS=/nologo /O2 /Qopenmp /Qiopenmp
    set MOD_FLAG=/module:
    set OSUFFIX=.obj
)

REM ソースディレクトリ
set UNST_DIR=src

REM ソースファイルの拡張子
set FSUFFIX=.f90

REM ビルド用ディレクトリ (中間ファイルの出力先)
set BUILD_DIR=build

REM モジュールが含まれるファイル名(拡張子なし) - コンパイル順序に注意
set UNST_MODULE_FNAME=UNST_Mod UNST_Mod2 UNST_Write UNST_Read UNST_Read2 UNST_Sub UNST_Riv

REM メインプログラムが含まれるファイル名(拡張子なし)
set MAIN_FNAME=UNST_Main

REM 除外対象ファイル名
REM set EXCLUDE_FNAME=

REM 実行ファイル名
set TARGET=unst.exe

REM -------------------------------------------
REM  STAGE-0 : ディレクトリ準備
REM -------------------------------------------
echo STAGE-0 Directory Setup
set OBJ_DIR=%BUILD_DIR%\obj
set MOD_DIR=%BUILD_DIR%\mod

REM ディレクトリの確認・生成
if exist "%BUILD_DIR%" rmdir /s /q "%BUILD_DIR%"
mkdir "%OBJ_DIR%"
mkdir "%MOD_DIR%"

REM -------------------------------------------
REM  STAGE-1 : UNSTモジュールのコンパイル
REM -------------------------------------------
echo STAGE-1 Compiling UNST Modules

for %%f in (%UNST_MODULE_FNAME%) do (
    echo [Compiling Module] %UNST_DIR%\%%f%FSUFFIX%
    if /I "%COMPILER%"=="gfortran" (
        "%FC%" %OPTIONS% -c %MOD_FLAG%"%MOD_DIR%" -o "%OBJ_DIR%\%%f%OSUFFIX%" "%UNST_DIR%\%%f%FSUFFIX%"
    ) else (
        "%FC%" %OPTIONS% /c %MOD_FLAG%"%MOD_DIR%" /object:"%OBJ_DIR%\%%f%OSUFFIX%" "%UNST_DIR%\%%f%FSUFFIX%"
    )
    if errorlevel 1 (
        echo *** ERROR Compile Module: %UNST_DIR%\%%f%FSUFFIX% ***
        exit /b 1
    )
)

REM -------------------------------------------
REM  STAGE-2 : UNSTソースのコンパイル
REM -------------------------------------------
echo STAGE-2 Compiling UNST Sources

for %%f in (%UNST_DIR%\*%FSUFFIX%) do (
    set skip=false
    set fname=%%~nf
    
    REM モジュールファイルはスキップ
    for %%m in (%UNST_MODULE_FNAME%) do (
        if /I "!fname!"=="%%m" set skip=true
    )
    
    if "!skip!"=="false" (
        echo [Compiling] %%f
        if /I "%COMPILER%"=="gfortran" (
            "%FC%" %OPTIONS% -c %MOD_FLAG%"%MOD_DIR%" -I"%MOD_DIR%" -o "%OBJ_DIR%\!fname!%OSUFFIX%" "%%f"
        ) else (
            "%FC%" %OPTIONS% /c %MOD_FLAG%"%MOD_DIR%" /I"%MOD_DIR%" /object:"%OBJ_DIR%\!fname!%OSUFFIX%" "%%f"
        )
        if errorlevel 1 (
            echo *** ERROR Compile: %%f ***
            exit /b 1
        )
    )
)

REM -------------------------------------------
REM  STAGE-3 : メインプログラムのコンパイル
REM -------------------------------------------
echo STAGE-3 Compiling Main Program

echo [Compiling Main] %UNST_DIR%\%MAIN_FNAME%%FSUFFIX%
if /I "%COMPILER%"=="gfortran" (
    "%FC%" %OPTIONS% -c %MOD_FLAG%"%MOD_DIR%" -I"%MOD_DIR%" -o "%OBJ_DIR%\%MAIN_FNAME%%OSUFFIX%" "%UNST_DIR%\%MAIN_FNAME%%FSUFFIX%"
) else (
    "%FC%" %OPTIONS% /c %MOD_FLAG%"%MOD_DIR%" /I"%MOD_DIR%" /object:"%OBJ_DIR%\%MAIN_FNAME%%OSUFFIX%" "%UNST_DIR%\%MAIN_FNAME%%FSUFFIX%"
)
if errorlevel 1 (
    echo *** ERROR Compile Main: %MAIN_FNAME%%FSUFFIX% ***
    exit /b 1
)

REM -------------------------------------------
REM  STAGE-4 : リンク処理
REM -------------------------------------------
echo STAGE-4 Linking

set OBJ_LIST=
if /I "%COMPILER%"=="gfortran" (
    for %%f in (%OBJ_DIR%\*.o) do set OBJ_LIST=!OBJ_LIST! "%%f"
) else (
    for %%f in (%OBJ_DIR%\*.obj) do set OBJ_LIST=!OBJ_LIST! "%%f"
)

if /I "%COMPILER%"=="gfortran" (
    "%FC%" %OPTIONS% !OBJ_LIST! -o "%BUILD_DIR%\%TARGET%"
) else (
    "%FC%" %OPTIONS% !OBJ_LIST! /Fe:"%BUILD_DIR%\%TARGET%"
)
if errorlevel 1 (
    echo *** ERROR Link ***
    exit /b 1
)

echo ==== success build ====
echo .exe: %BUILD_DIR%\%TARGET%

pause