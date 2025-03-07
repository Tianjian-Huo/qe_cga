# MyProject

## Overview
This repository contains scripts and configurations to set up and optimize molecular structures using `setup.py`. The process involves copying necessary files from the `cmds` directory and modifying configurations based on user input.

## Usage
Run the `setup.py` script and provide a molecular formula:
```bash
python setup.py
```
For example, entering `Ca1Ba3` will:
- Copy `optimize.input`, `cga.pbs`, `cmd_cmds.py`, and `config.ini` from the `cmds` directory.
- Modify `config.ini` according to the input molecular formula.
- Create `qein.json` in the `search` directory.

## Created Example Structure
```
Ca1Ba3_example/
│── cga.pbs/     
│── cmds_cmd.py    
│── config.ini             # Useless but esential file
│── search/
│   ├── config.ini       # Configuration file for the optimization process
│   ├── qein.json        # Quantum ESPRESSO input parameters settings
│   ├── optimize.input        
```

## Explanation of `search/config.ini`
The `config.ini` file contains parameters for optimization. Below is a breakdown of its contents:

| Parameter         | Description                      |
|------------------|--------------------------------|
| `core_num`       | Number of CPU cores to use    |
| `qe_dir`         | Path to Quantum ESPRESSO software |
| `qe_pseudo_dir`  | Path to pseudopotential files |

## Explanation of `search/qein.json`
The `qein.json` file defines parameters for Quantum ESPRESSO calculations. Here are the key settings:
```json
{
    "forc_conv_thr": 0.001,
    "nstep": 30,
    "ecutwfc": 40.0,
    "ecutrho": 240.0,
    "electron_maxstep": 30,
    "conv_thr": "1e-6",
    "degauss": 0.02,
    "smearing": "gaussian",
    "mixing_beta": 0.5,
    "CELL_PARAMETERS": "offset",
    "offset": 10.0
}
```

### Suggested Parameter Adjustments
| Parameter         | Suggested Values |
|------------------|----------------|
| `conv_thr`       | `1e-6` or `1e-7` |
| `degauss`        | `0.02` or `0.001` |
| `smearing`       | `gaussian` or `mp` |
| `mixing_beta`    | `0.5` or `0.3` |
| `CELL_PARAMETERS` | `offset` or `10` |

## Contribution
Feel free to fork this repository and submit pull requests for improvements.

For any questions or collaboration inquiries, contact me at **leonehuo@gmail.com**.

