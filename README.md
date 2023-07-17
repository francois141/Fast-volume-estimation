# Polytope Volume Estimation, but faster

## How to benchmark
```sh
make run_optimization1
cd benchmark_infrastructure
# conda activate base
python3 benchmark_generator.py
```

- Every optimization is in its own folder inside `./optimizations/`.
- The `./benchmark_infrastructure/main.cpp` is the one that does the cycle counting. It includes the files `polytope_optimized.h` and `linalg_optimized.h`.
- Thus, to include a specific optimization, you have to make sure the Make rule specifies the right include directory: `-I ./optimizitation/optmization-xx-yy`

## Build

```sh
make # compiles to be fast
make debug # compiles with debug and fsanitize flags
```

## GLPK

To install on your system, follow the steps below:

```bash
cd /tmp
wget http://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz
tar -xzf glpk-5.0.tar.gz
cd glpk-5.0
./configure
make --jobs=4
make check
sudo make install
```

