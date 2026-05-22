<p align="center">
  <h1>RCG-REG</h1>
  <p align="center">Image registration accuracy and efficiency with advanced Riemannian optimization.</p>
  <p align="center">
    <img alt="Build Status" src="https://img.shields.io/badge/Build-Passing-brightgreen">
    <img alt="License" src="https://img.shields.io/badge/License-MIT-blue">
    <img alt="PRs Welcome" src="https://img.shields.io/badge/PRs-Welcome-brightgreen">
    <img alt="GitHub Stars" src="https://img.shields.io/github/stars/user/repo?style=social">
  </p>
</p>

---

## The Strategic "Why"

> Traditional image registration methods often struggle with robustly aligning images under complex deformations, leading to suboptimal accuracy, slow convergence, and sensitivity to initial conditions. This can severely impact applications in medical imaging, computer vision, and scientific analysis, where precise alignment is paramount.

RCG-REG leverages the power of Riemannian nonmonotone conjugate gradient algorithms to provide a significantly more robust, accurate, and efficient solution for 2D and 3D image registration. By operating on the intrinsic geometry of the problem space, it overcomes the limitations of Euclidean methods, delivering superior alignment performance for even the most challenging datasets.

---

## Key Features

*   🚀 **Enhanced Accuracy**: Achieve superior registration precision by utilizing optimization on Riemannian manifolds.
*   ⚡ **Optimized Performance**: Experience faster convergence rates and reduced computational overhead compared to traditional methods.
*   ✨ **Robust 2D & 3D Support**: Seamlessly register both two-dimensional and three-dimensional images with a unified, powerful framework.
*   ⚙️ **Nonmonotone Strategy**: Benefit from a more forgiving optimization approach that can escape local minima more effectively, improving global convergence.
*   🌐 **Geometric Awareness**: Exploit the intrinsic structure of image transformation spaces, leading to more natural and accurate deformations.

---

## Technical Architecture

RCG-REG is built upon a robust foundation designed for high-performance scientific computing and advanced mathematical optimization.

| Technology                  | Purpose                                                                | Key Benefit                                                                  |
| :-------------------------- | :--------------------------------------------------------------------- | :--------------------------------------------------------------------------- |
| **MATLAB**                  | Primary Development Language & Runtime                                 | Robust environment for numerical computing, algorithm prototyping, and visualization. |
| **Riemannian Geometry**     | Mathematical Foundation for Optimization                               | Enables optimization on curved manifolds, enhancing registration accuracy and robustness. |
| **Conjugate Gradient Methods** | Core Optimization Algorithm                                            | Efficiently finds optimal solutions for large-scale problems with fewer iterations. |
| **Nonmonotone Strategy**    | Optimization Enhancement for Global Convergence                        | Improves global convergence by allowing temporary increases in the objective function, helping to escape local minima. |

### Directory Structure

```
RCG-REG/
├── 2d/
│   └── (MATLAB scripts and functions for 2D registration)
├── 3d/
│   └── (MATLAB scripts and functions for 3D registration)
├── LICENSE
└── README.md
```

---

## Operational Setup

### Prerequisites

Ensure you have the following software installed and configured:

*   **MATLAB**: R2018a or newer.
*   **MATLAB Image Processing Toolbox**: Essential for image manipulation and core registration functionalities.
*   **MATLAB Optimization Toolbox**: Required for advanced optimization routines.

### Installation

To get RCG-REG up and running, follow these steps:

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/user/RCG-REG.git
    cd RCG-REG
    ```

2.  **Add to MATLAB Path**:
    Open MATLAB and execute the following commands in the command window to add RCG-REG and its subdirectories to your MATLAB path. This ensures all necessary functions are accessible.
    ```matlab
    addpath(genpath(pwd)); % Adds current directory and all subdirectories
    savepath;              % Saves the updated path for future MATLAB sessions
    ```

---

