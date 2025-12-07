
# SM_Correlation – Smart Meter Correlation Analysis  
Final B.Sc. Project – Electrical & Electronics Engineering  

## 📌 Overview
This MATLAB project implements a correlation-based analysis algorithm designed to infer **customer-to-transformer connectivity** in low-voltage (LV) distribution grids using smart meter time-series data.  

The goal is to determine which customers are connected to the same transformer by analyzing the similarity (correlation) between their voltage measurements over time.  

This work is part of a final B.Sc. project submitted to the Afeka College of Engineering.

---

## 📂 Project Structure
```
SM_Correlation/
│
├──project_root/            # Main root – runs the full analysis pipeline
│  main.m
│  setup_paths.m
│
├─config/
│    get_data_config.m
│
├─algorithms/                # Our algorithems 
│    run_energy_balance.m
│    run_optimization.m
│    run_ecpc.m
│
├─utils/                     # Plots, reports and output files
│    load_or_generate_data.m
│    preprocess_data.m
│    evaluate_results.m
│
└─data/                       # Input datasets (Smart Meter voltage samples)
└── README.md    
```

---

## 🧠 Methodology Summary
The algorithms works in several stages:

1. **Data Loading**  
2. **Signal Preprocessing**  
3. **Algorithms Matrix Computation**  
4. **Transformer Grouping**  
5. **Visualization**

---

## 🚀 How to Run
```matlab
addpath(genpath('SM_Correlation'));
main
```

---

## 🛠 Requirements
- MATLAB R2021a or newer  
- Statistics and Machine Learning Toolbox (recommended)

---

## 👥 Authors
**Yariv Shossberger & Omri Itzhaki**  
Electrical & Electronics Engineering  
Afeka College of Engineering  


---

## 📬 Contact
yarivshossberger@gmail.com
