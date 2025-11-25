# 🚀 Windows HPC: Run H2O + SSDM Overnight

## Quick Start for Windows HPC in VS Code

---

## ✅ STEP 1: Verify R is Available

```bash
# In Git Bash or PowerShell
R --version
```

Should show R version 4.5.1 or similar. If not, R needs to be in your PATH.

---

## ✅ STEP 2: Navigate to Project

```bash
cd ~/BTEH/BTEH
pwd  # Should show: /c/Users/HarinPc/BTEH/BTEH or similar
```

---

## ✅ STEP 3: Verify Environment

```bash
# Start R
R

# In R console:
renv::status()
# Should show: "No issues found"

# Test packages
library(terra)
library(h2o)
library(SSDM)

# Exit R
quit(save = "no")
```

---

## ✅ STEP 4: Create Overnight Run Script

Create a file called `run_overnight.bat`:

```batch
@echo off
REM BTEH Overnight Run Script for Windows HPC
REM This will run H2O and SSDM training for all species

echo ========================================
echo BTEH Overnight Pipeline
echo Start time: %date% %time%
echo ========================================

cd /d "%~dp0BTEH\scripts_RemoteSensing for E and C_"

echo.
echo === H2O Training Run A ===
Rscript 03_h2o_train.R --run A --mode REPRO
if %errorlevel% neq 0 (
    echo ERROR in H2O Run A
    pause
    exit /b %errorlevel%
)

echo.
echo === H2O Training Run B ===
Rscript 03_h2o_train.R --run B --mode REPRO
if %errorlevel% neq 0 (
    echo ERROR in H2O Run B
    pause
    exit /b %errorlevel%
)

echo.
echo === SSDM Training Run A ===
Rscript 04_ssdm_train.R --run A --mode REPRO
if %errorlevel% neq 0 (
    echo ERROR in SSDM Run A
    pause
    exit /b %errorlevel%
)

echo.
echo === SSDM Training Run B ===
Rscript 04_ssdm_train.R --run B --mode REPRO
if %errorlevel% neq 0 (
    echo ERROR in SSDM Run B
    pause
    exit /b %errorlevel%
)

echo.
echo ========================================
echo Pipeline Complete!
echo End time: %date% %time%
echo ========================================
echo.
echo Results saved in: results/
pause
```

---

## ✅ STEP 5: Run Overnight

### **Option A: Double-click the .bat file**
1. Save `run_overnight.bat` in `~/BTEH/` folder
2. Double-click it before you leave
3. Window will stay open showing progress
4. Leave computer running overnight

### **Option B: Run from VS Code Terminal**

```bash
# In VS Code terminal (PowerShell or Git Bash)
cd ~/BTEH
./run_overnight.bat
```

### **Option C: Run in background (PowerShell)**

```powershell
# Start in background
cd ~/BTEH
Start-Process -FilePath "run_overnight.bat" -WindowStyle Minimized

# Or redirect output to log file
Start-Process -FilePath "run_overnight.bat" -RedirectStandardOutput "pipeline.log" -RedirectStandardError "errors.log" -NoNewWindow
```

---

## ✅ STEP 6: Monitor Progress

### **Check Log Files**

```bash
# View H2O log
tail -f BTEH/logs/03_h2o_train_A.log

# View SSDM log
tail -f BTEH/logs/04_ssdm_train_A.log
```

### **Check Results Folder**

```bash
# See what's been generated
ls -R BTEH/results/H2O/
ls -R BTEH/results/SSDM/
```

### **In VS Code**
- Open `BTEH/logs/` folder in Explorer
- Watch log files update in real-time
- Refresh to see new result files appear

---

## ✅ STEP 7: Alternative - Run Each Script Separately

If you want more control:

```bash
cd ~/BTEH/BTEH/"scripts_RemoteSensing for E and C_"

# Run one at a time
Rscript 03_h2o_train.R --run A --mode REPRO
# Wait for completion (2-4 hours)

Rscript 03_h2o_train.R --run B --mode REPRO
# Wait for completion (2-4 hours)

Rscript 04_ssdm_train.R --run A --mode REPRO
# Wait for completion (3-6 hours)

Rscript 04_ssdm_train.R --run B --mode REPRO
# Wait for completion (3-6 hours)
```

---

## ⏱️ Timeline for Tonight

| Task | Duration | Cumulative |
|------|----------|------------|
| H2O Run A (3 species) | 2-3 hours | 2-3 hours |
| H2O Run B (6 species) | 3-4 hours | 5-7 hours |
| SSDM Run A (3 species) | 2-3 hours | 7-10 hours |
| SSDM Run B (6 species) | 3-4 hours | 10-14 hours |

**Total: 10-14 hours** (perfect for overnight!)

---

## 🎯 What You'll Have Tomorrow Morning

```
results/
├── H2O/
│   ├── A/
│   │   ├── E3A/ ✅
│   │   ├── E4A/ ✅
│   │   └── E5A/ ✅
│   └── B/
│       ├── E1B/ ✅
│       ├── E2B/ ✅
│       ├── E3B/ ✅
│       ├── E4B/ ✅
│       ├── E5B/ ✅
│       └── E6B/ ✅
│
└── SSDM/
    ├── A/
    │   ├── E3A/ ✅
    │   ├── E4A/ ✅
    │   └── E5A/ ✅
    └── B/
        ├── E1B/ ✅
        ├── E2B/ ✅
        ├── E3B/ ✅
        ├── E4B/ ✅
        ├── E5B/ ✅
        └── E6B/ ✅
```

**All 9 species × 2 methods = 18 complete model sets!** 🎉

---

## 🐛 Troubleshooting

### **Issue: "Rscript not found"**
```bash
# Find R installation
where R
where Rscript

# Add to PATH if needed (PowerShell)
$env:Path += ";C:\Program Files\R\R-4.5.1\bin"
```

### **Issue: Script stops at first error**
- Check log files in `BTEH/logs/`
- Look for error messages
- Fix and restart from that step

### **Issue: Computer goes to sleep**
```powershell
# Prevent sleep (run before starting)
powercfg /change standby-timeout-ac 0
powercfg /change monitor-timeout-ac 30
```

---

## 💡 Pro Tips for Overnight Run

1. **Disable sleep mode** - Keep computer awake
2. **Close other programs** - Free up RAM
3. **Check disk space** - Need ~3-5 GB free
4. **Start early evening** - Give it 12-14 hours
5. **Check once before bed** - Make sure first script started
6. **Don't close terminal** - Keep window open

---

## ✅ Pre-Flight Checklist

Before you start the overnight run:

- [ ] R is installed and in PATH
- [ ] `renv::status()` shows no issues
- [ ] All packages load successfully
- [ ] Test run completed successfully
- [ ] Disk space > 5 GB free
- [ ] Computer won't sleep overnight
- [ ] `run_overnight.bat` created and saved
- [ ] Current directory is `~/BTEH`

---

## 🚀 START COMMAND (Copy-Paste This)

```bash
# Navigate to project
cd ~/BTEH

# Start overnight run
./run_overnight.bat
```

**Then leave it running overnight!** 🌙

---

## 📊 Check Progress in the Morning

```bash
# See what completed
ls -R BTEH/results/

# Check logs for any errors
grep -i error BTEH/logs/*.log

# Count completed models
find BTEH/results/H2O -name "prediction_*.tif" | wc -l
# Should show: 9

find BTEH/results/SSDM -name "ESDM_*.tif" | wc -l
# Should show: 9
```

---

**Good luck with your overnight run! You'll have all models ready by morning!** 🎯🐘
