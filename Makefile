CC=gcc
CPP=g++
M1_CC=gcc-12

# EXTRA_LIB_FLAGS = ##-static -zmuldefs

BOOST_INCLUDE_FLAG = -I /usr/include

##Lapack for Armadillo to do Cholesky decomposition
##GLPK for linear programming procedure
LIB_FLAGS = -larmadillo -lglpk -llapack -lcblas -lm -lgfortran $(EXTRA_LIB_FLAGS)
# LIB_FLAGS = -llapack -lblas -lglpk -lgfortran $(EXTRA_LIB_FLAGS)

OPT = -O3 -g -ffast-math -mfma # prevent the undefined reference error. not needed in debug bcs its included in -fsanitize
OPT_O2 = -O2 -mfma # O3 and ffast-math are not always the best
OPT_M1 = -O3 -ffast-math -march=native

SIMD = -march=native -mavx

DEBUG = -O0 -g -fsanitize=address,undefined
LINALG_FLAGS = -lglpk -llapack -lcblas -lm
ARMADILLO_FLAGS = -lstdc++ -larmadillo

# FINAL = -DARMA_NO_DEBUG
CCFLAGS_FAST = $(BOOST_INCLUDE_FLAG) $(OPT)
CCFLAGS_FAST_O2 = $(BOOST_INCLUDE_FLAG) $(OPT_O2)
CCFLAGS_FAST_M1 = -I /opt/homebrew/Cellar/armadillo/12.0.1/include $(FINAL) $(OPT_M1)
CCFLAGS_DEBUG = $(BOOST_INCLUDE_FLAG) $(DEBUG)

fast:
	$(CC) $(CCFLAGS_FAST) -I ./optimizations/optimization-fastest -o ./build/PolyVestFast main.c $(LINALG_FLAGS)
fast-m1:
	$(M1_CC) $(CCFLAGS_FAST_M1) -o ./build/PolyVestFast main.c polytope.c
debug:
	$(CC) $(CCFLAGS_DEBUG) -I ./optimizations/optimization-fastest -o ./build/PolyVestDebug main.c $(LINALG_FLAGS)
test_m1:
	$(M1_CC) $(CCFLAGS_FAST_M1) -L/opt/homebrew/Cellar/armadillo/12.0.1/lib -larmadillo -o ./build/PolyVestFastTest ./tests/test.cpp polytope.c
test_walk:
	$(CPP) $(CCFLAGS_DEBUG) -I ./optimizations/optimization-baseline -o ./build/PolyVestTestWalk ./tests/test_walk.cpp  ./tests/original_version/vol.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
test_matrix:
	$(CPP) $(CCFLAGS_DEBUG) -I ./optimizations/optimization-baseline -o ./build/PolyVestTestMatrix ./tests/test_matrix_operations.cpp  $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
test_cholesky:
	$(CPP) $(CCFLAGS_DEBUG) -I ./optimizations/optimization-baseline -o ./build/PolyVestTestCholesky ./tests/test_cholesky.cpp  $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
test_preprocess:
	$(CPP) $(CCFLAGS_DEBUG) -I ./optimizations/optimization-baseline -o ./build/PolyVestTestPreprocess ./tests/test_preprocess.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
generate_plots_original:
	$(CPP) $(OPT) -I ./optimizations/optimization-1-header-only -DARMA_NO_DEBUG -o ./benchmark_infrastructure/PolyVestGeneratePlots ./benchmark_infrastructure/main_orig.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp  $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
test_volume:
	$(CPP) $(CCFLAGS_DEBUG) -I ./optimizations/optimization-baseline -o ./build/PolyVestTestVolume ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
run_all_tests:
	make test_walk && ./build/PolyVestTestWalk
	make test_matrix && ./build/PolyVestTestMatrix
	make test_cholesky && ./build/PolyVestTestCholesky
	make test_preprocess && ./build/PolyVestTestPreprocess
	make test_volume && ./build/PolyVestTestVolume
run_optimization_baseline:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-baseline -o ./build/PolyVestOptimizationBaseline ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimizationBaseline baseline
run_optimization_fastest:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-fastest -o ./build/PolyVestOptimizationFastest ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimizationFastest fastest
run_optimization1:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-1-header-only -o ./build/PolyVestOptimization1 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization1 ours_o3
run_optimization2:
	$(CPP) $(CCFLAGS_FAST_O2) -I ./optimizations/optimization-2-o2 -o ./build/PolyVestOptimization2 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization2 ours_o2
run_optimization3:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-3-precompute-uball -o ./build/PolyVestOptimization3 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization3 ours_o4
run_optimization4:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-4-explicit-matrix-vector-mult -o ./build/PolyVestOptimization4 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization4 explicit-matrix-vector-mult
run_optimization5:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-5-scalar-replacement -o ./build/PolyVestOptimization5 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization5 scalar-replacement
run_optimization6:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-6-xoshiro-random -o ./build/PolyVestOptimization6 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization6 xoshiro-random
run_optimization7:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-7-remove-Adimensions -o ./build/PolyVestOptimization7 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization7 no_Adimensions
run_optimization8:
	$(CPP) $(CCFLAGS_FAST) -I ./optimizations/optimization-8-update-just-bounds -o ./build/PolyVestOptimization8 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization8 precompute_bounds
run_optimization9:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-9-simd -o ./build/PolyVestOptimization9 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization9 simd-9
run_optimization10:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-10-unrolling -o ./build/PolyVestOptimization10 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization10 unrolling-10
run_optimization11:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-11-simd-estimate -o ./build/PolyVestOptimization11 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization11 simd-11-estimate
run_optimization12:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-12-simd-estimate-padded -o ./build/PolyVestOptimization12 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization12 simd-12-estimate-padded
run_optimization13:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -march=skylake-avx512 -I ./optimizations/optimization-13-avx-512 -o ./build/PolyVestOptimization13 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization13 simd-13
run_optimization14:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-14-replace-div-with-mult -o ./build/PolyVestOptimization14 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization14 simd-14-replace-div-with-mult
run_optimization15:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-15-extract-mask -o ./build/PolyVestOptimization15 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization15 simd-15-extract-mask
run_optimization17:
	$(CPP) $(CCFLAGS_FAST) $(SIMD) -I ./optimizations/optimization-17-unrolling -o ./build/PolyVestOptimization17 ./benchmark_infrastructure/main.cpp $(LINALG_FLAGS)
	./build/PolyVestOptimization17 simd-17-unrolling

test_optimization_11:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-11-simd-estimate -o ./build/PolyVestTestOptimization11 ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization11
test_optimization_12:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-12-simd-estimate-padded -o ./build/PolyVestTestOptimization12 ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization12
test_optimization_13:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-13-avx-512  -o ./build/PolyVestTestOptimization13 ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization13
test_optimization_14:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-14-replace-div-with-mult  -o ./build/PolyVestTestOptimization14 ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization14
test_optimization_15:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-15-extract-mask  -o ./build/PolyVestTestOptimization15 ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization15

test_optimization_16:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-16-alignment-fix -o ./build/PolyVestTestOptimization16 ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization16

test_optimization_17:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-17-unrolling -o ./build/PolyVestTestOptimization17 ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization17

test_optimization_fastest:
	$(CPP) $(CCFLAGS_DEBUG) $(SIMD) -I ./optimizations/optimization-fastest -o ./build/PolyVestTestOptimization ./tests/test_volume.cpp ./tests/original_version/vol.cpp ./tests/original_version/walk.cpp $(LINALG_FLAGS) $(ARMADILLO_FLAGS)
	./build/PolyVestTestOptimization

build_for_execute_baseline:
	$(CC) $(CCFLAGS_FAST) -I ./optimizations/optimization-baseline -o ./build/PolyVestBaseline main.c $(LINALG_FLAGS)
build_for_execute_fastest:
	$(CC) $(CCFLAGS_FAST) -I ./optimizations/optimization-fastest -o ./build/PolyVestFastest main.c $(LINALG_FLAGS)
build_for_execute_basic_opt:
	$(CC) $(CCFLAGS_FAST) -I ./optimizations/optimization-6-xoshiro-random -o ./build/PolyVestBasicOpt main.c $(LINALG_FLAGS)
build_for_execute_algo_opt:
	$(CC) $(CCFLAGS_FAST) -I ./optimizations/optimization-8-update-just-bounds -o ./build/PolyVestAlgoOpt main.c $(LINALG_FLAGS)



clean:
	rm -rf ./build/*