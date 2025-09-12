# FireX: A Research Branch of Fire Dynamics Simulator (FDS)

FireX is a research branch of the Fire Dynamics Simulator (FDS) focused on integrating high-performance computing (HPC) functionalities. We try to keep FireX synchronized with the master branch of FDS, ensuring that any new commits to the master branch are incorporated into FireX within a few days.

### Key Features of FireX:

1. **GPU-accelerated Pressure Poisson Solver**
2. **VTK Output in HDF5 Format** – In addition to Smokeview output, enabling visualization in ParaView or other VTK-compatible tools.

---

## Compiling FireX

The compilation process for FireX is similar to FDS, as described in the [FDS Wiki](https://github.com/firemodels/fds/wiki). Below are the steps:

### **Step 1: Setup Your Environment**

Follow the instructions to set up your system environment:
  
  - [Setting up Windows environment](https://github.com/firemodels/fds/wiki/Setting-up-Windows-Environment)
  - [Setting up macOS environment](https://github.com/firemodels/fds/wiki/Setting-up-macOS-Environment)
  - [Setting up Linux environment](https://github.com/firemodels/fds/wiki/Setting-up-Linux-environment)

Note: CMake version 3.21 or higher is required for FireX.

#### **Additional Environment Variables for GPU Compilation**

Depending on your GPU type, set the following environment variables:

- **NVIDIA (CUDA):**
  ```sh
  export CUDA_DIR=<Path to CUDA installation (contains lib64)>
  export CUDA_MATH_DIR=<Path to CUDA math libraries (contains lib64)>
  ```
- **AMD (HIP):**
  ```sh
  export ROCM_DIR=<Path to ROCM installation (contains lib64)>
  ```
- **Intel (SYCL):**
  ```sh
  export SYCL_DIR=<Path to oneapi installation (contains lib)>
  ```

### **Step 2: Clone the Required Repositories**

Pick a local directory where you will clone `fds`, `hypre`, `sundials`, and 'hdf5'.  Optionally, you can name this directory `$FIREMODELS` in your startup script (e.g., `~/.bashrc`).  For example, 

```sh
export FIREMODELS=~/firemodels
```

Clone the repositories:

```sh
cd $FIREMODELS
git clone git@github.com:firemodels/fds.git
git clone git@github.com:hypre-space/hypre.git
git clone git@github.com:LLNL/sundials.git
git clone git@github.com:HDFGroup/hdf5.git
```

### **Step 3: Build FireX**

For compilation using **GNU compilers** use the GNU OpenMPI Linux target:
```sh
cd fds/Build/ompi_gnu_linux
git checkout FireX
./make_fds.sh --with-gpu=<cuda | hip>
```

For compilation using **Intel compilers** use the Intel Linux target:
```sh
cd fds/Build/impi_intel_linux
./make_fds.sh --with-gpu=sycl
```

---

## Verification Cases

- **GPU Verification:** `fds/Verification/test_gpu.fds`
- **VTK-HDF5 Verification:** `fds/Verification/VTK/beam_detector.fds`

---

## Example Environment Setup for Supercomputing Clusters

### **[Vista-TACC](https://tacc.utexas.edu/systems/vista/)** (NVIDIA - Grace Hopper)

To compile FireX using **GNU compilers**, add the following to your `~/.bash_profile`:

```sh
module load gcc/13.2.0
module load openmpi/5.0.7
module load cuda
module load nvidia_math
export FIREMODELS_FC=mpif90
export MPICH_DIR=$MPI_ROOT
export CUDA_DIR=$TACC_CUDA_DIR
export CUDA_MATH_DIR=$TACC_NVIDIA_MATH_DIR
export HYPRE_ENABLE_GPU_AWARE_MPI=ON
export NO_M64_FLAG=ON
```

### [Stampede-TACC](https://tacc.utexas.edu/systems/stampede3/) (INTEL Data Center GPU Max 1550s)
To compile FireX using INTEL compilers, add the following to your ~/.bash_profile:
```sh
module load intel/24.0
export MPICH_DIR=$MPI_ROOT
export HYPRE_ENABLE_GPU_AWARE_MPI=ON
export SYCL_DIR=/opt/intel/oneapi/compiler/2024.0
```
### [Aurora-ALCF](https://www.anl.gov/aurora) (Intel Data Center GPU Max Series)
To compile FireX using INTEL compilers, add the following to your ~/.bash_profile:
```sh
module load cmake
export FIREMODELS_CC=icx
export FIREMODELS_CXX=dpcpp
export FIREMODELS_FC=mpifort
export HYPRE_ENABLE_GPU_AWARE_MPI=ON
export SYCL_DIR=/opt/aurora/24.180.3/updates/oneapi/compiler/eng-20240629
```

### **[Polaris-ALCF](https://www.alcf.anl.gov/polaris)** (NVIDIA - A100)

To compile FireX using **GNU compilers**, add the following to your `~/.bash_profile`:

```sh
module use /soft/modulefiles
module load spack-pe-base cmake
module load PrgEnv-gnu
module load nvhpc-mixed
module load craype-accel-nvidia80
export FIREMODELS_CC=cc
export FIREMODELS_CXX=CC
export FIREMODELS_FC=ftn
export CUDA_DIR=$NVIDIA_PATH/cuda/12.2
export CUDA_MATH_DIR=$NVIDIA_PATH/math_libs/12.2
export HYPRE_ENABLE_GPU_AWARE_MPI=ON
```

### **[Frontier-OLCF](https://www.olcf.ornl.gov/frontier/)** (AMD - MI250)

To compile FireX using **GNU compilers**, add the following to your `~/.bash_profile`:

```sh
module load cmake/3.27.9
module load PrgEnv-gnu
module load rocm
module load craype-accel-amd-gfx90a
export FIREMODELS_CC=cc
export FIREMODELS_CXX=CC
export FIREMODELS_FC=ftn
export ROCM_DIR=$ROCM_PATH
export HYPRE_ENABLE_GPU_AWARE_MPI=ON
```

---

## Contact & Contributions

If you would like to contribute to FireX or report issues, feel free to open a GitHub issue or reach out to the maintainers via the FDS repository.






