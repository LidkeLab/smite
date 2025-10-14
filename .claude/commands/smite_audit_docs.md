# /smite_audit_docs

Audit LLM documentation for common issues, API errors, and testability problems.

## Purpose

Systematically check the smite LLM documentation for:
- **Untested code examples** that may contain bugs
- **API usage errors** (incorrect function calls, missing parameters)
- **Context-dependent requirements** not clearly documented
- **Inconsistent patterns** across documentation
- **Missing prerequisite setup** (e.g., DataROI for gaussBlobImage)

This command helps identify issues that can creep into LLM-generated documentation, especially patterns that get replicated across multiple files.

## Audit Categories

### 1. Code Example Validation

Check all code blocks in documentation:
- Extract MATLAB code snippets
- Identify incomplete setups (missing variable initialization)
- Find function calls with potentially missing prerequisites
- Flag examples that haven't been tested

**Common issues:**
- Missing `SMF.Data.DataROI` before `gaussBlobImage`
- Undefined variables referenced in code
- Function calls with wrong argument count/types
- Paths or file references that don't exist

### 2. API Consistency Checks

Verify function usage matches source code:
- Function signatures match implementation
- Required vs optional parameters correctly documented
- Default values are accurate
- Return values properly described

**Common issues:**
- "Optional" parameters that are context-dependent
- Default values incorrectly transcribed
- Outdated API from old examples

### 3. Context-Dependent Requirements

Find fields/parameters that are:
- Optional in some contexts
- Required in others
- Not clearly marked as context-dependent

**Example:**
- `SMF.Data.DataROI` is optional for analysis workflows
- But REQUIRED for `gaussBlobImage` image generation
- Should be documented as "Required for image generation"

### 4. Pattern Replication Issues

Find patterns that appear multiple times:
- Check if all instances are correct
- Identify copy-paste propagation of bugs
- Flag inconsistent usage of same API

**Common issues:**
- Bug in one example replicated to multiple files
- Inconsistent parameter names across examples
- Different patterns for the same operation

### 5. Missing Prerequisites

Code examples that jump into operations without setup:
- GPU initialization
- Camera calibration files
- Directory structures
- Required SMF fields

## Audit Process

### Step 1: Scan Documentation

```matlab
% Extract all code blocks from documentation
doc_files = find all *.md files in doc/llm-guide/
for each file:
    extract code blocks marked with ```matlab
    store with source file and line number
```

### Step 2: Static Analysis

For each code example:
- **Variable tracking**: Identify all variable uses and definitions
- **Function calls**: Extract function name and arguments
- **Prerequisites**: Check for required setup before operations
- **Dependencies**: Track cross-file dependencies

### Step 3: API Verification

Cross-reference with source code:
- Load actual function definitions
- Compare signatures
- Check default values
- Verify return values

### Step 4: Pattern Analysis

Group similar code patterns:
- Find all uses of specific functions
- Compare setup code
- Flag inconsistencies
- Suggest standardization

### Step 5: Generate Report

Create audit report with:
- **Summary**: Total issues by category
- **Critical issues**: Bugs that will prevent code from running
- **Warnings**: Potentially problematic patterns
- **Suggestions**: Improvements for clarity/consistency

## Report Format

```
LLM Documentation Audit Report
==============================
Generated: 2025-10-14

Summary
-------
Total files audited: 57
Total code blocks: 143
Issues found: 12

By Severity:
  CRITICAL: 2  (code won't run)
  MAJOR:    4  (incorrect results)
  MINOR:    6  (clarity/style)

By Category:
  Code Examples:           5
  API Consistency:         3
  Context-Dependent:       2
  Pattern Replication:     1
  Missing Prerequisites:   1

Critical Issues
---------------

1. workflows/spt-tracking.md:117
   SEVERITY: CRITICAL
   CATEGORY: Missing Prerequisites

   Issue: gaussBlobImage called without required DataROI

   Code:
   ```matlab
   SMF = smi_core.SingleMoleculeFitting();
   SMF.Fitting.PSFSigma = 1.3;
   [~, sequence] = smi_sim.GaussBlobs.gaussBlobImage(SPTSim.SMD, SMF);
   ```

   Problem: SMF.Data.DataROI must be set for gaussBlobImage

   Fix:
   ```matlab
   SMF = smi_core.SingleMoleculeFitting();
   SMF.Data.DataROI = [1, 1, 128, 128];  % Required
   SMF.Fitting.PSFSigma = 1.3;
   [~, sequence] = smi_sim.GaussBlobs.gaussBlobImage(SPTSim.SMD, SMF);
   ```

   Files to check:
   - examples/tracking-diffusion.md (CORRECT - has DataROI)
   - examples/multi-channel.md (CORRECT - has DataROI)
   - workflows/spt-tracking.md (BUG - missing DataROI)

2. core-concepts/smf-structure.md:210
   SEVERITY: CRITICAL
   CATEGORY: API Consistency

   Issue: PSFSigma default value incorrect

   Documentation says: 1.0
   Source code says: 1 (integer)

   Impact: Type mismatch, though functionally equivalent

   Fix: Change doc to match source (1 not 1.0)

Major Issues
------------

3. how-to/simulate-data.md:87
   SEVERITY: MAJOR
   CATEGORY: Code Examples

   Issue: Undefined variable 'SZ' used

   Code:
   ```matlab
   for nn = 1:NFrames
       N = poissrnd(Rho * SZ * SZ);  % SZ not defined
   ```

   Fix: Add SZ definition before loop

[... continue for all issues ...]

Pattern Analysis
----------------

gaussBlobImage Usage Pattern:
  Total occurrences: 7
  Correct setup: 6/7 (85.7%)
  Missing DataROI: 1/7 (14.3%)

  Recommendation: Add note to gaussBlobImage docs emphasizing
  DataROI requirement

Context-Dependent Requirements
-------------------------------

Fields marked "optional" but context-dependent:
1. SMF.Data.DataROI
   - Optional: For file-based analysis (loaded from data)
   - Required: For image generation (gaussBlobImage)
   - Recommendation: Document as "Required for image generation"

2. SMF.Data.CameraGain
   - Optional: For simulated data (default = 1)
   - Required: For real camera data (must match calibration)
   - Recommendation: Document as "Required for experimental data"

Recommendations
---------------

High Priority:
1. Fix 2 CRITICAL bugs immediately
2. Test all 143 code blocks
3. Add DataROI requirement to gaussBlobImage docs
4. Create pre-submission code validation script

Medium Priority:
5. Standardize gaussBlobImage usage pattern across all docs
6. Add "Context" section to SMF field docs (when required/optional)
7. Create troubleshooting doc for "optional but not really" fields

Low Priority:
8. Add inline comments to complex examples
9. Create cross-reference from context-dependent fields to examples
10. Consider automated testing of doc code blocks

Automated Fixes Available
--------------------------

The following issues can be automatically fixed:
- Issue #1: Add SMF.Data.DataROI line
- Issue #2: Change 1.0 to 1
- Issue #3: Add SZ definition

Run: /smite_fix_doc_issues

Manual Review Needed
--------------------

The following issues require human review:
- API design consistency (Issue #8)
- Terminology standardization (Issue #9)
- Example complexity vs clarity trade-offs

Test Coverage
-------------

Recommendation: Create automated test suite for doc examples

Proposed: MATLAB/tests/test_documentation_examples.m
- Extract code from each doc file
- Run in isolated environment
- Verify no errors
- Check outputs match expectations
- Run as part of CI/CD

Estimated effort: 4-6 hours to implement
Maintenance: ~15 min per doc update

Conclusion
----------

Overall documentation quality: GOOD (95/100)

The LLM-generated documentation is comprehensive and mostly
accurate. The 12 issues found represent edge cases and subtle
requirements that are easy to miss even for human authors.

Key takeaway: Context-dependent requirements need explicit
documentation to prevent confusion.

Next steps:
1. Apply automated fixes
2. Manually review remaining issues
3. Update status file with corrections
4. Consider API improvements to reduce confusion
```

## Usage Instructions

Run the audit:
```bash
/smite_audit_docs
```

Audit specific category:
```bash
/smite_audit_docs --category code-examples
/smite_audit_docs --category api-consistency
/smite_audit_docs --category context-dependent
```

Audit specific files:
```bash
/smite_audit_docs --files workflows/*.md
/smite_audit_docs --files examples/tracking-diffusion.md
```

Generate detailed report:
```bash
/smite_audit_docs --detailed
```

Apply automated fixes:
```bash
/smite_audit_docs --fix
```

## Implementation Notes

The audit should:
1. **Parse all markdown** in doc/llm-guide/
2. **Extract code blocks** (```matlab...```)
3. **Analyze statically** (don't execute yet)
4. **Cross-reference** with source code
5. **Generate report** with actionable fixes
6. **Update status file** with audit timestamp

Future enhancements:
- Actually execute code blocks in sandboxed environment
- Use MATLAB Code Analyzer for deeper checks
- Generate test suite from examples
- Track audit history over time
- Integrate with CI/CD pipeline

## Related Commands

- `/smite_validate_llm_docs` - Run comprehensive validation (includes testing)
- `/smite_build_llm_docs` - Build new documentation
- `/smite_update_llm_docs` - Update docs after code changes
- `/smite_docs_status` - Check documentation status

## Exit Criteria

Audit is successful when:
- All documentation files scanned
- All code blocks analyzed
- Report generated with issues categorized
- Automated fixes identified
- Status file updated

Audit reveals problems when:
- CRITICAL issues found (code won't run)
- MAJOR issues exceed threshold (>5%)
- Patterns show systematic problems
- API inconsistencies discovered
