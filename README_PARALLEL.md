# FreeSASA Parallel

[![C Tests](https://img.shields.io/badge/C%20tests-54%2F54%20passed-brightgreen.svg)](#)
[![Python Tests](https://img.shields.io/badge/Python%20tests-24%2F24%20passed-brightgreen.svg)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C Standard](https://img.shields.io/badge/C%20standard-C99-blue.svg)](#)
[![OpenMP](https://img.shields.io/badge/OpenMP-enabled-blue.svg)](#)

**FreeSASA Parallel** is a high-performance fork of the excellent
[FreeSASA](https://github.com/mittinatten/freesasa) library by
**Simon Mitternacht**. It adds OpenMP multi-threading and AVX-512/AVX2
SIMD vectorization to accelerate solvent accessible surface area
calculations, with a focus on MD trajectory analysis.

The public C and Python APIs are **identical** to the original.
Existing code that uses FreeSASA will work without any changes.

---

## Why Parallelism Matters

The original FreeSASA is accurate, easy to use, and well-tested.
However, modern computational biology workflows regularly process
**molecular dynamics trajectories** with thousands of frames and
systems of tens of thousands of atoms.

- A single 58,000-atom frame takes ~1 second serially.
- A 10,000-frame trajectory would take ~3 hours without parallelism.
- With this fork, the same trajectory finishes in **under 30 minutes**
  on an 8-core workstation.

Two levels of parallelism are combined:

1. **Atom-level (OpenMP):** The inner S&R and L&R loops are distributed
   across CPU cores using `#pragma omp parallel for`.
2. **SIMD (AVX-512 / AVX2):** Neighbor distance checks in the S&R
   algorithm process 16 atoms simultaneously per CPU cycle using
   float32 Structure-of-Arrays (SoA) layout and AVX-512 intrinsics,
   with automatic fallback to AVX2 or scalar on older hardware.
3. **Frame-level:** The new `freesasa_calc_structures_parallel()` API
   (and the Python `calcStructuresParallel()` wrapper) processes multiple
   trajectory frames concurrently, giving near-linear speedup with core count.

---

## Performance

All benchmarks measured on a system with AVX-512 support.
Test structure: **1AON (GroEL-GroES complex, 58,674 atoms)**.

### Single-Structure Speedup

#### Lee-Richards (L&R) — standard precision (20 slices)

| Threads | Time    | Speedup |
| :-----: | :-----: | :-----: |
| 1       | 1.185 s | 1.00×   |
| 2       | 0.641 s | 1.85×   |
| 4       | 0.371 s | 3.19×   |
| 8       | 0.230 s | **5.15×** |
| 16      | 0.183 s | **6.48×** |

#### Lee-Richards (L&R) — high precision (100 slices)

| Threads | Time    | Speedup |
| :-----: | :-----: | :-----: |
| 1       | 4.078 s | 1.00×   |
| 4       | 1.116 s | 3.65×   |
| 8       | 0.618 s | **6.60×** |
| 16      | 0.483 s | **8.44×** |

#### Shrake-Rupley (S&R) — 100 test points, with AVX-512 SIMD

| Threads | Time    | Speedup |
| :-----: | :-----: | :-----: |
| 1       | 0.266 s | 1.00×   |
| 2       | 0.191 s | 1.39×   |
| 4       | 0.149 s | 1.79×   |
| 8       | 0.133 s | 2.00×   |
| 16      | 0.126 s | **2.11×** |

> S&R scaling is limited by the neighbor-list construction (O(N²)) which
> dominates at low point counts. At higher resolutions (500+ test points)
> the SIMD inner loop contributes more and scaling improves further.

---

### Trajectory Speedup (frame-level parallel API)

#### Large system — 8 frames × 58,674 atoms

| Mode                    | Time    | Speedup |
| :---------------------- | :-----: | :-----: |
| Serial (1 frame at a time) | 9.38 s | 1.00× |
| Parallel — 2 frames     | 4.84 s  | 1.94×   |
| Parallel — 4 frames     | 2.52 s  | 3.72×   |
| Parallel — 8 frames     | 1.38 s  | **6.79×** |

#### Typical MD system — 32 frames × 602 atoms (1UBQ)

| Mode                    | Time    | Speedup |
| :---------------------- | :-----: | :-----: |
| Serial (1 frame at a time) | 0.365 s | 1.00× |
| Parallel — 2 frames     | 0.168 s | 2.17×   |
| Parallel — 4 frames     | 0.087 s | 4.21×   |
| Parallel — 8 frames     | 0.052 s | 7.07×   |
| Parallel — 16 frames    | 0.044 s | **8.36×** |

---

## What Was Changed

All changes are backward-compatible. No existing API was modified or removed.

| Component           | Change |
| :------------------ | :----- |
| `src/sasa_lr.c`     | pthread → OpenMP; removed thread cap; dynamic scheduling |
| `src/sasa_sr.c`     | AVX-512/AVX2 SIMD inner loop; float32 SoA neighbor cache |
| `src/nb.c`          | Two-phase parallel neighbor-list construction (was serial O(N²)) |
| `src/freesasa.c`    | New `freesasa_calc_structures_parallel()` function |
| `src/freesasa.h`    | New `freesasa_calc_structures_parallel()` declaration |
| `src/main.cc`       | CLI auto-detects core count via `omp_get_max_threads()` |
| `CMakeLists.txt`    | CMake build system (replaces autotools for the parallel build) |
| `config.h`          | Compile-time feature flags; `#ifndef` guards for safe inclusion |

### Portability
- **C standard:** C99 throughout. `aligned_alloc` (C11) replaced with
  a `posix_memalign` / `_aligned_malloc` portable wrapper.
- **No-SIMD fallback:** When AVX-512/AVX2 are not detected at compile
  time, the code falls back to `#pragma omp simd` auto-vectorization.
- **Single-threaded mode:** `OMP_NUM_THREADS=1` or `parameters.n_threads=1`
  restores fully serial behavior with no conflicts.

---

## Quick Start

### Build

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
ctest --output-on-failure   # 54/54 tests should pass
```

### Control Thread Count

```bash
# Environment variable (affects all OpenMP programs)
export OMP_NUM_THREADS=8

# Or inline for a single run
OMP_NUM_THREADS=8 freesasa --n-threads=8 structure.pdb
```

### Python — Trajectory Analysis

```python
import freesasa

params = freesasa.Parameters()
params.setNThreads(8)  # 8 frames processed simultaneously

results = freesasa.calcStructuresParallel(frame_structures, params)
total_areas = [r.totalArea() for r in results]
```

---

## Acknowledgment

This fork is based entirely on the work of **Simon Mitternacht**.
All credit for the core algorithm, API design, and test suite belongs to him.

If you use this software in research, please cite the original paper:

> Mitternacht S. **FreeSASA: An open source C library for solvent accessible
> surface area calculations.** *F1000Research.* 2016;5:189.
> doi:[10.12688/f1000research.7931.1](https://doi.org/10.12688/f1000research.7931.1)
