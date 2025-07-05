# makefile
# coded by k.yamamura
#

# コンパイラの指定
FC = gfortran

# コンパイルオプションの指定
FFLAGS = -g -O2 -cpp -ffree-line-length-none -fopenmp
# debugモード
# FFLAGS = -g -O0 -cpp -ffree-line-length-none -fopenmp -fbacktrace -fcheck=all

# リンクオプションの指定
LDFLAGS = -fopenmp

# ソースファイルの拡張子
FSUFFIX = .f90

# ソースファイルの取得
SOURCES = $(wildcard *$(FSUFFIX))

# オブジェクトファイルの拡張子
OSUFFIX = .o

# オブジェクトファイルの取得
OBJECTS = $(patsubst %$(FSUFFIX),%$(OSUFFIX),$(SOURCES))

# 実行ファイル名
TARGET = unst.exe

# デフォルトターゲット
all: clean $(TARGET)

# 実行ファイルの生成
$(TARGET): $(OBJECTS)
	$(FC) $(LDFLAGS) -o $@ $^ || (rm -f $(OBJECTS) $(SOURCES:.f90=.d); exit 1)

# オブジェクトファイルの生成
%$(OSUFFIX): %$(FSUFFIX)
	$(FC) $(FFLAGS) -c $< || (rm -f $@ $(patsubst %$(FSUFFIX),%.d,$<); exit 1)

# 依存関係ファイルの生成
%.d: %$(FSUFFIX)
	$(FC) -cpp -MM $< > $@ || (rm -f $@; exit 1)

# 依存関係ファイルのインクルード
-include $(SOURCES:.f90=.d)

# クリーンアップ
.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TARGET) $(SOURCES:.f90=.d) $(SOURCES:.f90=.mod)

# 実行後の後処理
.PHONY: post_process
post_process:
	rm -f $(OBJECTS) $(SOURCES:.f90=.d) $(SOURCES:.f90=.mod)

# デフォルトターゲットの実行後に後処理を実行
all: $(TARGET) post_process