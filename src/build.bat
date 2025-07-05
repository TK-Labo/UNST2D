@echo off
setlocal enabledelayedexpansion

REM コンパイラ設定
set FC=gfortran

REM フラグ設定
set FCFLAGS_RELEASE=-O2 -g -cpp -ffree-line-length-none -fopenmp
set FCFLAGS_DEBUG=-O0 -g -cpp -Wall -Wextra -fcheck=all -fbacktrace -ffree-line-length-none -fopenmp

REM デフォルトのモード
set FCFLAGS=%FCFLAGS_RELEASE%

REM 実行ファイル名
set EXEC=unst.exe

REM 自動クリーンアップするかどうか（0=しない、1=する）
set AUTO_CLEAN=1

REM コマンドライン引数のチェック
if "%1"=="debug" (
    echo build by mode debug
    set FCFLAGS=%FCFLAGS_DEBUG%
) else if "%1"=="release" (
    echo build by mode release
    set FCFLAGS=%FCFLAGS_RELEASE%
) else if "%1"=="clean" (
    echo ready clean up
    del /Q *.o %EXEC% *.mod 2>NUL
    exit /b 0
) else if "%1"=="no-clean" (
    echo build without auto cleanup
    set AUTO_CLEAN=0
) else (
    if not "%1"=="" (
        echo unknown option: %1
        echo use method: build.bat [debug^|release^|clean^|no-clean]
        exit /b 1
    )
    echo default - build by mode release
)

REM f90ファイルのリストを取得
set "SOURCES="
for %%F in (*.f90) do (
    set "SOURCES=!SOURCES! %%F"
)

REM ソースファイルがない場合のチェック
if "!SOURCES!"=="" (
    echo error: NO found .f90 files.
    exit /b 1
)

REM 既存のオブジェクトファイルをクリーンアップ
del /Q *.o %EXEC% *.mod 2>NUL

REM インクルードファイルがあるか確認して処理
for %%F in (*.f90) do (
    findstr /i /c:"include " %%F > nul
    if not errorlevel 1 (
        echo Found include directives in: %%F
    )
)

REM ファイルをコンパイル
echo ready start compile...
for %%F in (*.f90) do (
    findstr /i /c:"module " %%F > nul
    if not errorlevel 1 (
        echo compile module: %%F
        %FC% %FCFLAGS% -c %%F
        if errorlevel 1 (
            echo compile error: %%F
            exit /b 1
        )
    )
)

REM 残りのファイルをコンパイル
echo next compling...
set "OBJECTS="
for %%F in (*.f90) do (
    findstr /i /c:"module " %%F > nul
    if errorlevel 1 (
        echo compiling: %%F
        %FC% %FCFLAGS% -c %%F
        if errorlevel 1 (
            echo compile error: %%F
            exit /b 1
        )
    )
    set "OBJECTS=!OBJECTS! %%~nF.o"
)

REM 実行ファイルの生成
echo ready start link
%FC% %FCFLAGS% -o %EXEC% %OBJECTS%
if errorlevel 1 (
    echo link error
    exit /b 1
)

echo success build: %EXEC%

REM ビルド後のクリーンアップ（オブジェクトファイルを削除）
if %AUTO_CLEAN% EQU 1 (
    echo cleaning object files...
    del /Q *.o *.mod 2>NUL
    echo cleanup complete, keeping only executable file.
)

exit /b 0