# Context-Dependent Requirements in smite

## The Problem

Many SMF (SingleMoleculeFitting) fields are documented as "optional" with default values (often `[]`), but are actually **required in specific contexts**. This leads to:

1. **Confusing errors** - Missing required field causes cryptic downstream error
2. **Documentation ambiguity** - "Optional" doesn't tell full story
3. **User frustration** - Code works in one context, fails in another
4. **Testing gaps** - Field never tested in context where it's required

## Examples of Context-Dependent Fields

### 1. SMF.Data.DataROI

**Documented as:**
```matlab
SMF.Data.DataROI = [];  % Optional, default []
```

**Reality:**
- **Optional** when loading from files (gets size from data)
- **REQUIRED** when calling `gaussBlobImage()` (needs output dimensions)
- **Optional** for most analysis workflows
- **Required** for image generation workflows

**Better documentation:**
```matlab
SMF.Data.DataROI = [YStart, XStart, YEnd, XEnd];
% [1, 1, 128, 128] for full 128x128 image
%
% REQUIRED for:
%   - smi_sim.GaussBlobs.gaussBlobImage() - defines output size
%   - Image generation from SMD
%
% OPTIONAL for:
%   - File-based analysis (loaded from data)
%   - LocalizeData workflows (gets from image stack)
%
% If omitted when required: Causes index error in image generation
```

**Detection in code:**
```matlab
% WRONG - will fail
SMF = smi_core.SingleMoleculeFitting();
[~, img] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);  % ERROR!

% RIGHT - includes required context
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.DataROI = [1, 1, 128, 128];  % Required for this operation
[~, img] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);  % Works!
```

### 2. SMF.Data.CameraGain

**Documented as:**
```matlab
SMF.Data.CameraGain = 1;  % Default 1
```

**Reality:**
- **Acceptable default** for simulated data (unit gain)
- **Must match calibration** for real camera data
- **Critical** for photon number accuracy
- **Less critical** for position-only analysis

**Better documentation:**
```matlab
SMF.Data.CameraGain = 1;  % electrons/ADU
%
% Simulated data: Use 1 (default)
% Real sCMOS: Load from calibration file (typical: 0.5-5)
% Real EMCCD: Set to EM gain value (typical: 50-300)
%
% Impact if incorrect:
%   - Photon numbers wrong (affects SNR, precision estimates)
%   - Positions usually OK (gain cancels in ratio)
%   - Thresholding may fail (wrong photon-based thresholds)
```

### 3. SMF.Data.PixelSize

**Documented as:**
```matlab
SMF.Data.PixelSize = 0.1;  % micrometers, default 0.1
```

**Reality:**
- **Any value works** for pixel-space analysis
- **Must be accurate** for physical units
- **Required for calibration** (camera µm to physical µm)
- **Affects clustering** if distance-based

**Better documentation:**
```matlab
SMF.Data.PixelSize = 0.1;  % micrometers per pixel
%
% Common values:
%   - 100nm pixel, 100x objective, no magnifier: 0.1 µm
%   - 16µm pixel, 100x, 1.6x magnifier: 0.1 µm
%   - 6.5µm pixel, 60x: 0.108 µm
%
% Used for:
%   - Converting pixel coordinates to physical units
%   - Drift correction distance thresholds
%   - Clustering distance parameters
%   - Physical diffusion coefficients
%
% If incorrect:
%   - Wrong physical scale (everything scaled incorrectly)
%   - Drift correction may fail (wrong search radius)
%   - Clustering may under/over-connect
```

## General Pattern

### When a field is context-dependent:

1. **Don't just say "optional"** - Specify contexts
2. **List when required** - Explicit use cases
3. **List when optional** - Where defaults work
4. **Document failure mode** - What happens if missing/wrong
5. **Provide examples** - Show correct usage in both contexts

### Documentation Template

```matlab
SMF.Category.FieldName = <default_value>;  % <units>
%
% REQUIRED for:
%   - <operation 1> - <why needed>
%   - <operation 2> - <why needed>
%
% OPTIONAL for:
%   - <operation A> - <why default works>
%   - <operation B> - <why not needed>
%
% DEFAULTS to:
%   <default_value> - <what this means>
%
% EXAMPLES:
%   <common value 1> - <use case>
%   <common value 2> - <use case>
%
% IMPACT if incorrect:
%   - <consequence 1>
%   - <consequence 2>
%
% SEE ALSO: <related fields>, <related docs>
```

## API Design Implications

### Current Issues

1. **Silent failures** - Missing field causes downstream error
2. **Inconsistent validation** - Some functions check, others don't
3. **Poor error messages** - Index errors instead of "missing DataROI"
4. **No hints** - Function doesn't tell you what's needed

### Improvements to Consider

#### Option 1: Validate at function entry

```matlab
function [Model, Data] = gaussBlobImage(SMD, SMF, Bg, Density)
    % Validate required fields
    if isempty(SMF.Data.DataROI)
        error(['gaussBlobImage requires SMF.Data.DataROI to be set.\n' ...
               'Example: SMF.Data.DataROI = [1, 1, 128, 128];']);
    end

    % ... rest of function
end
```

**Pros:** Clear error message, fail fast
**Cons:** Adds validation code to every function

#### Option 2: Separate setup functions

```matlab
% Instead of:
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.DataROI = [1, 1, 128, 128];
SMF.Fitting.PSFSigma = 1.3;

% Provide:
SMF = smi_core.SingleMoleculeFitting.forImageGeneration(128, 128);
% Automatically sets DataROI and other required fields
```

**Pros:** Clearer intent, groups related settings
**Cons:** More API surface, duplicate code

#### Option 3: Context-aware validation

```matlab
function SMF = validateForContext(SMF, context)
    % context = 'analysis' | 'simulation' | 'image_generation'

    switch context
        case 'image_generation'
            assert(~isempty(SMF.Data.DataROI), ...
                'DataROI required for image generation');
        case 'real_data'
            assert(~isempty(SMF.Data.CalibrationFile), ...
                'Calibration required for real data');
        % ... etc
    end
end
```

**Pros:** Explicit context checking
**Cons:** User must remember to call, adds complexity

#### Option 4: Documentation improvements only

Keep API as-is, but improve documentation:
- Mark fields as "Context-Dependent"
- Add "Required For" sections
- Improve error messages where possible
- Add validation examples to docs

**Pros:** No API changes, backward compatible
**Cons:** Still relies on user reading docs

## Recommended Approach

**Short term** (for v1.3.0 documentation):
1. Add "Context-Dependent Requirements" section to docs
2. Mark affected fields explicitly in SMF structure docs
3. Add validation examples to function docs
4. Improve error messages where possible (non-breaking)

**Medium term** (v1.4.0):
5. Add input validation to high-risk functions
6. Create helper constructors for common contexts
7. Add warnings for suspicious configurations
8. Create troubleshooting guide for context errors

**Long term** (v2.0.0):
9. Consider context-aware SMF variants
10. Redesign API to make requirements explicit
11. Add comprehensive validation framework
12. Breaking changes OK in major version

## Audit Checklist

When reviewing documentation for context-dependent requirements:

- [ ] Field documented with "optional" or "default"?
- [ ] Check: Is it truly optional in ALL contexts?
- [ ] List contexts where field is required
- [ ] List contexts where field is optional
- [ ] Document what happens if missing when required
- [ ] Add example code showing correct setup
- [ ] Check if function validates input (add if not)
- [ ] Add error message improvement if needed
- [ ] Cross-reference other docs using same field
- [ ] Test code examples in both contexts

## Testing Strategy

Create tests that verify context-dependent behavior:

```matlab
function test_gaussBlobImage_requires_DataROI()
    % This should FAIL with helpful error
    SMF = smi_core.SingleMoleculeFitting();
    % Deliberately omit DataROI

    SMD = smi_core.SingleMoleculeData.createSMD();
    SMD.X = [64]; SMD.Y = [64];
    SMD.Photons = [1000]; SMD.Bg = [5];
    SMD.PSFSigma = [1.3]; SMD.FrameNum = [1];

    try
        [~, ~] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);
        error('Should have failed without DataROI');
    catch ME
        % Check error message is helpful
        assert(contains(ME.message, 'DataROI'), ...
            'Error should mention DataROI requirement');
    end
end

function test_gaussBlobImage_works_with_DataROI()
    % This should SUCCEED
    SMF = smi_core.SingleMoleculeFitting();
    SMF.Data.DataROI = [1, 1, 128, 128];  % Required!

    SMD = smi_core.SingleMoleculeData.createSMD();
    SMD.X = [64]; SMD.Y = [64];
    SMD.Photons = [1000]; SMD.Bg = [5];
    SMD.PSFSigma = [1.3]; SMD.FrameNum = [1];

    [Model, Data] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);

    assert(~isempty(Model), 'Should return model');
    assert(size(Model, 1) == 128, 'Size from DataROI');
    assert(size(Model, 2) == 128, 'Size from DataROI');
end
```

## Summary

**Key Principles:**

1. **Be explicit** - Don't hide requirements in context
2. **Document both** - When required AND when optional
3. **Show examples** - Correct usage in both contexts
4. **Fail helpfully** - Clear error messages
5. **Test both** - Required context AND optional context

**For LLM Documentation:**

When generating examples:
- Always check if field is context-dependent
- Include required setup for the example context
- Add comments explaining why field is needed
- Cross-reference to full explanation
- Test code actually runs in stated context

**For Code Review:**

Red flags:
- Function assumes field is set without checking
- Error message doesn't mention missing field
- Documentation says "optional" without context
- Example code missing required setup
- No tests for missing required field

## Related Documents

- `.claude/commands/smite_audit_docs.md` - How to audit for these issues
- `.claude/commands/smite_validate_llm_docs.md` - Comprehensive validation
- `doc/llm-guide/core-concepts/smf-structure.md` - SMF field documentation
- `doc/llm-guide/troubleshooting/common-mistakes.md` - User-facing errors

## Version History

- 2025-10-14: Initial documentation of pattern
- Future: Update as API changes or more patterns discovered
