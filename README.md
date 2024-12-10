# Lactation Curve Parameter Estimation

This repository contains Python code to estimate and analyze lactation curve parameters for dairy cows. The code implements two methods to evaluate and compare lactation parameters based on Wood's lactation model, integrating adjustments for various factors such as region, year, parity, and milking frequency. The accompanying paper details the methodology and results *(provide DOI once online)*.

## Overview

### Key Features
1. **Wood's Lactation Model**: Implements Wood's curve to estimate daily milk yield over the lactation period.
2. **Two Methods**: 
   - **Method 1**: Base lactation curve parameters with adjustments from a predefined dictionary based on region, year, and milking frequency.
   - **Method 2**: Fine-tuned lactation curve parameters using annual herd-level milk production (AHMP) and parity-specific structures.
3. **Graphical Output**: Plots lactation curves for different parities under each method.

## Requirements

- Python 3.8+
- Required Python packages:
  - `numpy`
  - `pandas`
  - `scipy`
  - `matplotlib`

Install the required packages using:
```bash
pip install numpy pandas scipy matplotlib
```

## Usage

### Functions

#### 1. `calculate_lactation_group_yield`
Calculates lactation group yield based on annual herd-level milk production (AHMP), the number of milking cows, and parity structure percentages.

#### 2. `get_t_values`
Generates an array of DIM (days in milk) values for plotting the lactation curve.

#### 3. `get_y_values_wood_curve`
Calculates daily milk yield based on Wood's lactation model parameters.

#### 4. `calc_integral_wood_curve`
Computes the integral of Wood's lactation curve over the lactation period (305 days).

#### 5. `get_wood_parameters`
Retrieves or fine-tunes parameters \(a, b, c\) for Wood's curve based on adjustments and input variables.

#### 6. `method1`
Generates lactation curves for different parities using pre-adjusted parameters. It corresponds to the Method1 in the paper. 

#### 7. `method2`
Generates lactation curves for different parities by fine-tuning parameters based on AHMP and parity structure percentages. It corresponds to the Method2 in the paper. 

### Running the Code

#### Method 1 Example:
```python
milking_frequency = '3x/d'
region = 'New York'
method1('2016', region, milking_frequency)
```

#### Method 2 Example:
```python
parity_structure_pct = (0.359525, 0.263485, 0.3769925)
num_milking_cows = 1284
AHMP = 18405332   # Unit: kg
method2('2016', 'New York', '3x/d', AHMP, parity_structure_pct, num_milking_cows)
```

### Outputs
- **Console Output**: Displays parameter values (\(a, b, c\)) and calculated MY\(_{305d}\) values.
- **Graphs**: Plots lactation curves for different parities.

## Adjustment Dictionary
The `adjustment_dict` provides adjustments for Wood's model parameters based on a study of Li et al. (2022):
- **Parity**: Parity 1, 2, or 3+.
- **Year**: 2006 to 2016.
- **Month**: January to December.
- **Region**: U.S. regions like Northeast, West Coast, etc.
- **Milking Frequency**: `2x/d` or `3x/d`.

## Citation
Li, M., G.J.M. Rosa, K.F. Reed, and V.E. Cabrera. 2022. Investigating the effect of temporal, geographic, and management factors on US Holstein lactation curve parameters. J Dairy Sci 105:7525–7538. doi:10.3168/jds.2022-21882.

---

For any questions, please contact [Yijing Gong: gong44@wisc.edu] or [Haowen Hu: hh598@cornell.edu].
