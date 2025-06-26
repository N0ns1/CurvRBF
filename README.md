# CurvRBF: Mean Curvature-Controllable RBFs

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the official C++ implementation for the paper:

**CurvRBF: Mean Curvature-Controllable Radial Basis Functions for Implicit Geological Modeling**
*Yuxiang Chen, Hao Deng, Yang Zheng, Wenwen Shi, Xiancheng Mao*

*(Note: The paper is currently under review. A link to the publication will be added upon acceptance.)*


---

## 🚀 Key Features

*   **Mean Curvature Control**: Directly manipulate the local shape of implicit surfaces by setting curvature constraints.
*   **Enhanced HRBF**: Extends the traditional Hermite Radial Basis Function (HRBF) method with novel mean curvature constraints.
*   **Flexible Constraints**: Apply constraints directly on the surface (`on-contact`) or in the surrounding field (`off-contact`) for fine-grained control.
*   **Improved Accuracy with Sparse Data**: Achieves better reconstruction results, especially when input data is limited.

---

## 🛠️ Prerequisites

Before you begin, ensure you have the following installed:

*   **IDE**: **Visual Studio 2019** (v16.10 or newer) or **Visual Studio 2022**.
*   **C++ Standard**: The project requires a compiler with C++20 support.
*   **CUDA Toolkit**: **(Optional)** Required *only* if you want to use the GPU-accelerated Dual Marching Cubes (DMC,https://github.com/rogrosso/tmc) surfacer. The core `CurvRBF` logic does not depend on CUDA.

## ⚙️ Building from Source

Follow these steps to compile the project:

1.  Clone the repository:
    ```bash
    git clone https://github.com/N0ns1/CurvRBF.git
    cd CurvRBF
    ```
2.  Open the `CurvRBF.sln` file with Visual Studio.
3.  In the top toolbar, set the solution configuration to **Release** and the platform to **x64**.
    
4.  Build the solution by selecting `Build > Build Solution` from the menu or by pressing `Ctrl+Shift+B`.
5.  The executable `CurvRBF.exe` will be generated in the `x64/Release` directory.

## ▶️ Usage

The program is a command-line tool. The basic syntax is:

```bash
CurvRBF.exe <input_file> [output_path]
```

*   `<input_file>`: Path to the input point cloud file (e.g., `.xyz` format with normals).
*   `[output_path]`: (Optional) The directory where output files will be saved.

#### Example 1: Specifying an output directory

This will read the input file and save all results (e.g., `.obj` model, info file) into `D:\output\folder`.

```bash
.\x64\Release\CurvRBF.exe C:\data\models\input.xyz D:\output\folder
```

#### Example 2: Omitting the output directory

If you don't provide an output path, the results will be saved in the same directory as the input file (`C:\data\models\`).

```bash
.\x64\Release\CurvRBF.exe C:\data\models\input.xyz
```

## 📄 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
