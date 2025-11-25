# IRL-Based-SDRE
**State-Dependent Riccati Equation Control using Integral Reinforcement Learning**

MATLAB implementation accompanying the paper  
**"Extension of the State-Dependent Riccati Equation Method using Integral Reinforcement Learning"**  
(Aria Rashidinejad Meibodi & Mahbod Gholamali Sinaki, 2024–2025)

This repository contains three clean, reproducible scripts that clearly compare:

| File              | Description                                                                                   | Key Feature                                             |
|-------------------|-----------------------------------------------------------------------------------------------|----------------------------------------------------------|
| `Uncontrolled.m`  | Open-loop simulation (u = 0)                                                                  | System is unstable → states blow up                      |
| `SDRE.m`          | Classical (online) SDRE                                                                       | Solves the state-dependent ARE **at every integration step** using `care` |
| `SDRE_IRL.m`      | Proposed IRL-based SDRE (our contribution)                                                    | Solves the ARE only once **offline** via Integral Reinforcement Learning (no `care` needed during simulation) |

### System
ẋ₁ = (1 − x₁²)x₁ + x₂ + u

ẋ₂ = (1 + x₁x₂)x₁ − x₂ + u

B = [1; 1], Q = I₂×₂, R = 1  
Initial condition x₀ = [3; 1], simulation 0–10 s

### Requirements
- MATLAB R2018b or newer
- Control System Toolbox → required only for `SDRE.m`

### How to run
```matlab
>> Uncontrolled   % explodes
>> SDRE           % converges (slow because it solves ARE thousands of times)
>> SDRE_IRL       % converges fast (ARE already learned offline)
