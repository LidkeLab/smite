# Documentation Fixes Needed

**Priority fixes** identified by automated testing audit on 2025-10-14.

---

## HIGH PRIORITY (Fix ASAP)

### 1. workflows/spt-tracking.md:117 - Missing DataROI

**Issue**: gaussBlobImage example missing required field

**Current code**:
```matlab
SMF = smi_core.SingleMoleculeFitting();
SMF.Fitting.PSFSigma = 1.3;
[~, sequence] = smi_sim.GaussBlobs.gaussBlobImage(SPTSim.SMD, SMF);
```

**Fixed code**:
```matlab
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.DataROI = [1, 1, 128, 128];  % Required for image generation
SMF.Fitting.PSFSigma = 1.3;
[~, sequence] = smi_sim.GaussBlobs.gaussBlobImage(SPTSim.SMD, SMF);
```

**Impact**: Users copying this code will get immediate errors.

---

### 2. api-reference/navigation.md:628 - Wrong method call pattern

**Issue**: Documentation shows static method call for instance method

**Current**:
```matlab
smi.SMLM.fullAnalysis()  % This doesn't work!
```

**Fixed**:
```matlab
% Create object first
SMLMobj = smi.SMLM(SMF);
% Then call instance method
SMLMobj.fullAnalysis()
```

**Impact**: Confuses users about static vs instance methods.

---

### 3. Multiple files - Case-sensitivity with SMF

**Issue**: `SMF` conflicts with `smf` from MATLAB fuzzy toolbox

**Locations**:
- api-reference/navigation.md:245
- api-reference/navigation.md:636
- api-reference/smi-core.md:911
- api-reference/smi-psf.md:92
- Others

**Current**:
```matlab
SMF = smi_core.SingleMoleculeFitting();
% ... later in doc ...
SMF.Data.DataSize = [128, 128];  % May fail if fuzzy toolbox loaded
```

**Fixed**:
```matlab
% Option 1: Use full path in examples
SMF_obj = smi_core.SingleMoleculeFitting();

% Option 2: Add note at top of files
% Note: SMF is case-sensitive. If you see errors about 'smf' from
% the fuzzy toolbox, use: clear SMF; SMF = smi_core.SingleMoleculeFitting();
```

**Impact**: Code fails for users with fuzzy toolbox loaded.

---

## MEDIUM PRIORITY (Fix Soon)

### 4. how-to/simulate-data.md:87 - Undefined variable

**Issue**: Example uses `SZ` without defining it

**Current**:
```matlab
for nn = 1:NFrames
    N = poissrnd(Rho * SZ * SZ);  % SZ not defined!
end
```

**Fixed**:
```matlab
SZ = SMF.Data.DataSize(1);  % Image size
for nn = 1:NFrames
    N = poissrnd(Rho * SZ * SZ);
end
```

---

### 5. reference/file-formats.md - Multiple encoding issues

**Issues**:
- Line 128: Invalid expression (encoding problem)
- Line 532: Invalid text character (non-ASCII)
- Line 543: Incorrect '=' operator usage

**Action**: Re-type these code examples from scratch (copy-paste may have introduced encoding issues)

---

### 6. Add context-dependent field notes

**Files to update**:
- core-concepts/smf-structure.md
- api-reference/smi-sim.md (gaussBlobImage section)
- how-to/simulate-data.md

**Add notes like**:
```markdown
### SMF.Data.DataROI

**Type**: [Y_min, X_min, Y_max, X_max] (1×4 vector)
**Default**: [] (empty, loaded from data file)
**Context**:
- **Optional** for analysis workflows (loaded from data)
- **REQUIRED** for image generation (gaussBlobImage, genRandomBlobImage)
- **REQUIRED** for simulation workflows

**Example**:
```matlab
% For simulation/image generation
SMF.Data.DataROI = [1, 1, 128, 128];  % Required
```

---

## LOW PRIORITY (Improve Clarity)

### 7. Standardize gaussBlobImage examples

**Pattern to use**:
```matlab
% Always show this complete setup
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.DataROI = [1, 1, 128, 128];  % Required for image generation
SMF.Data.DataSize = [128, 128];
SMF.Fitting.PSFSigma = 1.3;
[~, sequence] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);
```

**Files to check**:
- All examples using gaussBlobImage (7 total)
- Most are correct (6/7), ensure consistency

---

### 8. Add troubleshooting section

**Location**: Create `troubleshooting/context-dependent-fields.md`

**Content outline**:
```markdown
# Context-Dependent Fields

Some SMF fields are optional in one workflow but required in another.

## SMF.Data.DataROI

- Analysis: Optional (loaded from file)
- Simulation: **Required**
- Error message: "DataROI must be set"
- Solution: `SMF.Data.DataROI = [1, 1, SZ, SZ]`

## SMF.Data.CameraGain

- Simulation: Optional (default = 1)
- Real data: **Required** (from calibration)
- Error message: "Gain not specified"
- Solution: Load calibration file or set manually

[... more examples ...]
```

---

### 9. Fix classdef examples

**Issue**: Classdef blocks can't execute (need to be in .m files)

**Examples** (many in API reference):
- api-reference/smi-cluster.md:64
- api-reference/smi-core.md:139
- Many others

**Solution**: Either:
1. Mark as syntax-only (add comment: `% Class definition (syntax example)`)
2. Move to separate code fence type: ` ```matlab-syntax `
3. Add note: "This shows class structure. See source code for full implementation."

---

## Documentation Improvement Tasks

### Add to gaussBlobImage documentation

In `api-reference/smi-sim.md`, add prominent note:

```markdown
#### ⚠️ Important: DataROI Requirement

`gaussBlobImage` **requires** `SMF.Data.DataROI` to be set, even though
this field is optional for analysis workflows.

```matlab
% This will fail:
SMF = smi_core.SingleMoleculeFitting();
[~, img] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);  % Error!

% This works:
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.DataROI = [1, 1, 128, 128];  % Required!
[~, img] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);  % Success
```
```

---

### Add to SMF structure documentation

In `core-concepts/smf-structure.md`, add "Context" column to all field tables:

| Field | Type | Default | Context | Description |
|-------|------|---------|---------|-------------|
| `DataROI` | 1×4 | [] | **Required for simulation/image generation** | Region of interest |
| `CameraGain` | scalar | 1 | **Required for real data** | Camera gain (e/count) |
| ... | ... | ... | ... | ... |

---

## Testing Integration

### Add pre-commit hook

Create `.git/hooks/pre-commit`:

```bash
#!/bin/bash
# Test documentation examples before commit

echo "Testing documentation examples..."
cd MATLAB/tests
matlab -batch "test_doc_examples" > /dev/null 2>&1

if [ $? -ne 0 ]; then
    echo "Documentation tests failed! Please fix before committing."
    exit 1
fi

echo "Documentation tests passed!"
exit 0
```

---

## Automated Fix Script (Future)

Could create `MATLAB/tests/fix_common_doc_issues.m`:

```matlab
function fix_common_doc_issues()
% Automatically fix common documentation issues

    % Fix #1: Add DataROI to gaussBlobImage examples
    fix_missing_dataroi();

    % Fix #2: Replace SMF with explicit class path
    fix_smf_case_sensitivity();

    % Fix #3: Add context notes to field descriptions
    add_context_annotations();

    % Fix #4: Verify all examples against source code
    verify_api_examples();
end
```

---

## Summary Checklist

- [ ] Fix workflows/spt-tracking.md (DataROI)
- [ ] Fix api-reference/navigation.md (method call)
- [ ] Add SMF case-sensitivity warnings
- [ ] Fix how-to/simulate-data.md (SZ variable)
- [ ] Re-type file-formats.md examples (encoding)
- [ ] Add context-dependent field notes
- [ ] Standardize gaussBlobImage examples
- [ ] Create troubleshooting guide
- [ ] Update gaussBlobImage documentation
- [ ] Add Context column to SMF tables
- [ ] Set up pre-commit testing hook

**Estimated total time**: 4-6 hours

**Priority order**:
1. Fix critical bugs (1-2 hours)
2. Add context annotations (1-2 hours)
3. Standardize patterns (1-2 hours)
4. Set up testing infrastructure (1-2 hours)

---

*Generated by LLM Documentation Audit - 2025-10-14*
