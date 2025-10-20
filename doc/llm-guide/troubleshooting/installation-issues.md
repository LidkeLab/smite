---
title: "Installation Issues"
category: "troubleshooting"
level: "beginner"
tags: ["troubleshooting", "installation", "setup", "errors", "path", "toolbox"]
prerequisites: []
related: ["../getting-started/installation.md", "../getting-started/quickstart.md", "gpu-problems.md", "compilation-errors.md"]
summary: "Diagnose and resolve common installation and setup problems with smite"
estimated_time: "10-20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Installation Issues

## Purpose

This guide helps you diagnose and resolve common problems encountered during smite installation and initial setup. Each issue includes symptoms, diagnosis steps, and solutions.

## Prerequisites

- Basic familiarity with MATLAB
- Access to MATLAB R2021a or later
- Administrator/sudo access (for some solutions)

## Overview

Installation issues typically fall into these categories:

1. **PATH issues** - MATLAB can't find smite files
2. **MATLAB version incompatibility** - Version too old or incompatible features
3. **Missing toolboxes** - Required MATLAB toolboxes not installed
4. **Permission errors** - Insufficient access rights
5. **Startup script problems** - Configuration errors

This guide provides systematic approaches to identify and fix each type of problem.

## Issue 1: setupSMITE Not Found

### Symptoms

```matlab
>> setupSMITE
Undefined function or variable 'setupSMITE'.
```

Or when starting MATLAB:
```
Undefined function or variable 'setupSMITE'.
Error in startup (line X)
```

### Diagnosis

This means MATLAB can't find the setupSMITE.m file, indicating a PATH problem.

**Step 1: Check if setupSMITE exists**

```matlab
which setupSMITE
```

**Expected output if working:**
```
C:\Users\YourName\Documents\MATLAB\smite\MATLAB\setupSMITE.m
```

**If not found:**
```
'setupSMITE' not found.
```

**Step 2: Verify smite directory structure**

```matlab
% Check if directory exists
exist('C:\Users\YourName\Documents\MATLAB\smite\MATLAB', 'dir')
```

Should return `7` if directory exists.

**Step 3: Check current path**

```matlab
path
```

Look for smite\MATLAB in the list.

### Solutions

**Solution A: startup.m doesn't exist or is misconfigured**

1. Create or edit startup.m:

```matlab
% Open startup file (creates if doesn't exist)
edit(fullfile(userpath, 'startup.m'))
```

2. Add these lines (adjust path for your system):

**Linux/Mac:**
```matlab
addpath '~/Documents/MATLAB/smite/MATLAB'
setupSMITE
```

**Windows:**
```matlab
addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
setupSMITE
```

3. Save the file
4. Restart MATLAB

**Solution B: Wrong path in startup.m**

Check the path is correct:

```matlab
% View startup.m
type(fullfile(userpath, 'startup.m'))
```

Common mistakes:
- Typos in path (e.g., "smite\MATAB" instead of "smite\MATLAB")
- Missing quotes around paths with spaces
- Using forward slashes on Windows or vice versa
- Path points to wrong directory level

**Correct examples:**
```matlab
% Correct (Windows)
addpath 'C:\Users\John Smith\Documents\MATLAB\smite\MATLAB'

% Incorrect (missing quotes with spaces)
addpath C:\Users\John Smith\Documents\MATLAB\smite\MATLAB

% Correct (Linux/Mac)
addpath '~/Documents/MATLAB/smite/MATLAB'

% Incorrect (missing MATLAB subdirectory)
addpath '~/Documents/MATLAB/smite'
```

**Solution C: Manual path addition**

Temporary fix (until MATLAB restart):

```matlab
% Add path manually (adjust for your system)
addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
setupSMITE
```

**Solution D: Wrong smite location**

If smite is installed elsewhere, use the full path:

```matlab
% Find where smite actually is
[status, result] = system('where setupSMITE.m');  % Windows
% or
[status, result] = system('find ~ -name setupSMITE.m 2>/dev/null');  % Linux/Mac
```

Update startup.m with the correct path.

### Verification

After applying solution:

```matlab
which setupSMITE
```

Should show full path to setupSMITE.m

```matlab
setupSMITE
```

Should print:
```
SMITE version: 1.x.x
```

## Issue 2: MATLAB Version Incompatibility

### Symptoms

```matlab
Error using ...
Unrecognized function or variable 'arguments'.
```

Or:
```
This version of MATLAB is not supported.
Requires R2021a or later.
```

Or older MATLAB-specific errors like:
```
Error: File ... contains features that are not supported by MATLAB.
```

### Diagnosis

**Check MATLAB version:**

```matlab
version
```

**Example output:**
```
9.10.0.1684407 (R2021a)
```

**Check if version is too old:**

```matlab
verLessThan('matlab', '9.10')
```

Returns `1` if version is older than R2021a (too old for smite).

### Solutions

**Solution A: Update MATLAB**

If you have an active license:

1. Download latest MATLAB from: https://www.mathworks.com/downloads
2. Install update
3. Restart MATLAB
4. Verify: `version`

**Solution B: Use older smite release**

If you can't update MATLAB, try an older smite version:

1. Check smite releases: https://github.com/LidkeLab/smite/releases
2. Look for releases compatible with your MATLAB version (check release notes)
3. Download and install that version

**Note:** Older releases may lack features and bug fixes.

**Solution C: Upgrade MATLAB license**

If you have an expired or student license:

1. Contact your institution's IT department
2. Or purchase new license from MathWorks

### Verification

```matlab
version  % Should show R2021a or later
SMF = smi_core.SingleMoleculeFitting()  % Should work without errors
```

## Issue 3: Missing Toolboxes

### Symptoms

```matlab
Error using ...
Undefined function 'gpuArray' for input arguments of type 'double'.
```

Or:
```
Requires Image Processing Toolbox.
Error using imfilter
```

Or:
```
License checkout failed.
License Manager Error -4
Maximum number of users for Statistics_Toolbox reached.
```

### Diagnosis

**Step 1: Check installed toolboxes**

```matlab
ver
```

**Expected output includes:**
```
MATLAB                                                Version 9.10        (R2021a)
Image Processing Toolbox                              Version 11.3        (R2021a)
Parallel Computing Toolbox                            Version 7.4         (R2021a)
Statistics and Machine Learning Toolbox               Version 12.1        (R2021a)
```

**Step 2: Check specific toolbox**

```matlab
license('test', 'image_toolbox')
license('test', 'distrib_computing_toolbox')
license('test', 'statistics_toolbox')
```

Returns `1` if available, `0` if not.

**Step 3: Identify which toolbox is missing**

Common error patterns:

| Error Message | Missing Toolbox |
|--------------|-----------------|
| `gpuArray`, `gpuDevice` | Parallel Computing Toolbox |
| `imfilter`, `imadjust` | Image Processing Toolbox |
| `fitlm`, `kmeans`, `normrnd` | Statistics and Machine Learning Toolbox |
| `fit`, `cfit` | Curve Fitting Toolbox (optional) |
| `fmincon`, `lsqnonlin` | Optimization Toolbox (optional) |

### Solutions

**Solution A: Install missing toolboxes**

1. **Open Add-On Explorer:**
   - MATLAB Home tab → Add-Ons → Get Add-Ons

2. **Search for required toolbox:**
   - "Image Processing Toolbox"
   - "Parallel Computing Toolbox"
   - "Statistics and Machine Learning Toolbox"

3. **Install toolbox** (requires MathWorks account and license)

4. **Restart MATLAB**

**Solution B: License issues**

If toolbox is installed but license fails:

1. **Check license server** (for network licenses):
   ```matlab
   license('inuse')
   ```

2. **Contact your license administrator:**
   - May need to wait for license to free up
   - May need additional license seats

3. **Verify license file** (for standalone licenses):
   - Help → Licensing → Activate Software

**Solution C: Use alternative MATLAB installation**

If you have access to another MATLAB installation with required toolboxes:
- Use that installation for smite
- Or run smite on a university/company server with full toolboxes

### Verification

```matlab
% Verify required toolboxes
ver  % Should list all required toolboxes

% Test each toolbox
license('test', 'image_toolbox')
license('test', 'distrib_computing_toolbox')
license('test', 'statistics_toolbox')
% All should return 1

% Quick functional test
gpuDevice  % Tests Parallel Computing Toolbox (requires NVIDIA GPU)
```

## Issue 4: Permission Errors

### Symptoms

```
Error using save
Unable to write file ... : permission denied.
```

Or:
```
Error: Access is denied.
Cannot create directory: ...
```

Or when running tests:
```
Error using mkdir
Permission denied: C:\Program Files\MATLAB\...
```

### Diagnosis

**Step 1: Check where smite is installed**

```matlab
which setupSMITE
```

**Problem locations:**
- `C:\Program Files\` (Windows) - requires admin rights
- `/usr/local/` (Linux/Mac) - requires sudo
- System-protected directories

**Step 2: Check write permissions**

```matlab
% Try to create a file in smite directory
testFile = fullfile(fileparts(which('setupSMITE')), 'test.txt');
try
    fid = fopen(testFile, 'w');
    if fid == -1
        disp('Cannot write to smite directory')
    else
        fclose(fid);
        delete(testFile);
        disp('Write permissions OK')
    end
catch ME
    disp('Permission error:')
    disp(ME.message)
end
```

**Step 3: Check if running MATLAB as admin**

Windows: Check if MATLAB window title shows "Administrator"
Linux/Mac: Check if running as root (not recommended)

### Solutions

**Solution A: Move smite to user directory**

Best practice: Install in user-writable location.

1. **Move smite:**

**Windows:**
```
Move from: C:\Program Files\MATLAB\smite
To: C:\Users\YourName\Documents\MATLAB\smite
```

**Linux/Mac:**
```bash
mv /usr/local/smite ~/Documents/MATLAB/smite
```

2. **Update startup.m** with new path

3. **Restart MATLAB**

**Solution B: Fix directory permissions**

If you must keep smite in current location:

**Windows:**
1. Right-click on smite folder
2. Properties → Security → Edit
3. Add your user account
4. Grant "Full Control"
5. Apply to all subfolders

**Linux/Mac:**
```bash
# Make user owner of directory
sudo chown -R $USER ~/Documents/MATLAB/smite

# Set appropriate permissions
chmod -R u+w ~/Documents/MATLAB/smite
```

**Solution C: Run MATLAB as administrator** (Not recommended)

**Windows:**
- Right-click MATLAB icon
- "Run as administrator"

**Linux/Mac:**
```bash
sudo matlab
```

**Warning:** Running as admin/root can cause other permission issues. Only use temporarily.

**Solution D: Change output directories**

If permission errors occur during analysis:

```matlab
% Change where results are saved
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.ResultsDir = 'C:\Users\YourName\Documents\smite_results';
```

### Verification

```matlab
% Test write permissions
testDir = fullfile(fileparts(which('setupSMITE')), 'test_permissions');
[success, msg] = mkdir(testDir);
if success
    rmdir(testDir);
    disp('Permissions OK')
else
    disp(['Permission problem: ' msg])
end

% Run quick test
smi_core.LocalizeData.unitTest()  % Should complete without permission errors
```

## Issue 5: Startup Script Errors

### Symptoms

MATLAB starts but shows errors:

```
Error in startup (line X)
```

Or infinite loops/hangs during startup.

Or:
```
Error using run
Unable to read file 'startup.m'. No such file or directory.
```

### Diagnosis

**Step 1: Check startup.m syntax**

```matlab
edit(fullfile(userpath, 'startup.m'))
```

Look for:
- Missing quotes around paths
- Typos in function names
- Extra or missing semicolons (not usually critical but can cause issues)
- Invalid MATLAB syntax

**Step 2: Check what's in startup.m**

```matlab
type(fullfile(userpath, 'startup.m'))
```

**Step 3: Test startup.m manually**

```matlab
% Temporarily rename startup.m
movefile(fullfile(userpath, 'startup.m'), fullfile(userpath, 'startup_backup.m'))

% Restart MATLAB (should start without errors)

% Then run startup manually to see exact error
run(fullfile(userpath, 'startup_backup.m'))
```

### Solutions

**Solution A: Fix syntax errors**

Common mistakes and fixes:

**Missing quotes:**
```matlab
% Wrong
addpath C:\Users\John Smith\Documents\MATLAB\smite\MATLAB

% Correct
addpath 'C:\Users\John Smith\Documents\MATLAB\smite\MATLAB'
```

**Wrong slash direction:**
```matlab
% Wrong (Windows)
addpath 'C:/Users/John/Documents/MATLAB/smite/MATLAB'

% Correct (Windows)
addpath 'C:\Users\John\Documents\MATLAB\smite\MATLAB'
```

**Calling setup before adding path:**
```matlab
% Wrong order
setupSMITE
addpath 'C:\Users\John\Documents\MATLAB\smite\MATLAB'

% Correct order
addpath 'C:\Users\John\Documents\MATLAB\smite\MATLAB'
setupSMITE
```

**Solution B: Simplify startup.m**

Create minimal working version:

```matlab
% Minimal startup.m for smite
try
    addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
    setupSMITE
catch ME
    warning('smite setup failed: %s', ME.message)
end
```

The try-catch prevents startup errors from blocking MATLAB launch.

**Solution C: Reset startup.m**

If startup.m is corrupted:

1. **Backup current version:**
```matlab
copyfile(fullfile(userpath, 'startup.m'), fullfile(userpath, 'startup_old.m'))
```

2. **Create fresh startup.m:**
```matlab
edit(fullfile(userpath, 'startup.m'))
```

3. **Add only smite setup:**
```matlab
addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
setupSMITE
```

4. **Save and restart MATLAB**

**Solution D: Disable startup temporarily**

To diagnose if startup.m is the problem:

1. **Rename startup.m:**
```matlab
movefile(fullfile(userpath, 'startup.m'), fullfile(userpath, 'startup_disabled.m'))
```

2. **Restart MATLAB** (should start normally)

3. **Add smite manually:**
```matlab
addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
setupSMITE
```

4. **If smite works,** problem was in startup.m - recreate it correctly

### Verification

```matlab
% View final startup.m
type(fullfile(userpath, 'startup.m'))

% Should be simple and correct, e.g.:
% addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
% setupSMITE

% Restart MATLAB and verify
% Should see: SMITE version: 1.x.x

% Test basic functionality
SMF = smi_core.SingleMoleculeFitting()
```

## Additional Common Issues

### Issue: External Software Setup Fails

**Symptoms:**
```
Error in setupExternalSoftware
```

**Diagnosis:**
```matlab
which setupExternalSoftware
```

**Solution:**
Check that `ExternalSoftware` directory exists:

```matlab
exist(fullfile(fileparts(fileparts(which('setupSMITE'))), 'ExternalSoftware'), 'dir')
```

If missing, re-download/clone smite repository completely.

### Issue: Conflicting MATLAB Toolboxes

**Symptoms:**
```
Error: Name 'function_name' conflicts with existing function
```

**Diagnosis:**
```matlab
which function_name -all
```

Shows multiple locations for same function.

**Solution:**
```matlab
% Remove conflicting path
rmpath('/path/to/conflicting/toolbox')

% Make permanent by updating startup.m
```

### Issue: antivirus Blocking Files

**Symptoms:**
- Setup runs but functions don't work
- "File not found" errors for files that exist
- Intermittent failures

**Diagnosis:**
Check antivirus logs for blocked MATLAB/smite files

**Solution:**
1. Add smite directory to antivirus exceptions
2. Add MATLAB installation to antivirus exceptions
3. Temporarily disable antivirus during installation (not recommended long-term)

## Systematic Troubleshooting Workflow

If you're still having issues, follow this systematic approach:

### Step 1: Clean Install

```matlab
% 1. Remove smite path
rmpath(genpath('C:\Users\YourName\Documents\MATLAB\smite'))

% 2. Restart MATLAB

% 3. Delete and re-download smite
% (In file explorer or terminal, delete smite directory)
% Re-clone from GitHub

% 4. Add path and setup
addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
setupSMITE
```

### Step 2: Verify Prerequisites

```matlab
% Check MATLAB version
version  % Should be R2021a or later

% Check toolboxes
ver  % Should show Image Processing, Parallel Computing, Statistics

% Check GPU (if applicable)
gpuDevice  % Should show NVIDIA GPU info
```

### Step 3: Run Diagnostic Tests

```matlab
% Test 1: Basic functionality
SMF = smi_core.SingleMoleculeFitting()

% Test 2: Simple localization
B = smi_sim.GaussBlobs.genRandomBlobImage();
B = poissrnd(B);
LD = smi_core.LocalizeData(B, SMF);
[SMD] = LD.genLocalizations();

% Test 3: Unit test
smi_core.LocalizeData.unitTest()
```

### Step 4: Get Help

If still not working:

1. **Gather diagnostic information:**
```matlab
version
ver
gpuDevice  % if applicable
which setupSMITE
path
```

2. **Report issue:**
   - GitHub Issues: https://github.com/LidkeLab/smite/issues
   - Include MATLAB version, OS, error messages, diagnostic info
   - Include minimal reproducible example

## Prevention Tips

**1. Use recommended installation location:**
   - `~/Documents/MATLAB/smite` (all platforms)

**2. Keep startup.m simple:**
   - Only essential paths
   - Use try-catch for error handling

**3. Verify after updates:**
   - After MATLAB updates
   - After smite updates
   - After system updates

**4. Regular testing:**
```matlab
% Quick test after any changes
setupSMITE
SMF = smi_core.SingleMoleculeFitting()
```

**5. Document your configuration:**
   - Note any custom paths
   - Note any modified startup.m settings
   - Keep backups of working configurations

## See Also

- [Installation Guide](../getting-started/installation.md) - Complete setup instructions
- [Quick Start](../getting-started/quickstart.md) - Get running in 5 minutes
- [GPU Problems](gpu-problems.md) - GPU-specific troubleshooting
- [Compilation Errors](compilation-errors.md) - Mex/CUDA compilation issues

## Summary

Most installation issues fall into these categories:

1. **PATH issues** → Fix startup.m, verify paths
2. **Version problems** → Update MATLAB or use older smite
3. **Missing toolboxes** → Install required toolboxes
4. **Permissions** → Move to user directory or fix permissions
5. **Startup errors** → Simplify and fix startup.m syntax

The key to successful troubleshooting:
- Systematic diagnosis
- Isolate the problem
- Apply targeted solution
- Verify fix works

If problems persist after trying these solutions, report an issue on GitHub with complete diagnostic information.
