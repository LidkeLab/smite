# Validate LLM Documentation

Validate that LLM documentation is consistent, accurate, and complete.

## Purpose

Performs comprehensive validation checks on documentation:
- File existence and structure
- Cross-reference integrity
- Frontmatter consistency
- Link validation
- Code example syntax
- Optional: Code example execution
- Manifest accuracy

## Check Prerequisites

1. Verify `.claude/llm-docs-status.json` exists
   - If missing, suggest running `/smite_build_llm_docs` first
2. Verify `doc/llm-guide/` directory exists
3. Report what will be validated

## Validation Checks

### 1. File Existence

**Check:**
- All documents in manifest exist on disk
- All documents in status file exist on disk
- No orphaned files (exist but not in manifest)
- Required metadata files exist (manifest.json, index.md)

**Report:**
```
File Existence Check
====================
✓ All manifest documents exist (30/30)
✓ All status file documents exist (30/30)
✓ No orphaned documents found
✓ manifest.json exists
✓ index.md exists

Status: PASSED
```

**On Failure:**
```
File Existence Check
====================
✗ Missing documents (2):
  • workflows/roi-analysis.md (in manifest, not on disk)
  • examples/3d-localization.md (in status, not on disk)

⚠ Orphaned documents (1):
  • draft/experimental-feature.md (on disk, not in manifest)

✓ manifest.json exists
✓ index.md exists

Status: FAILED
Action: Remove missing entries from manifest/status or create missing files
```

### 2. Cross-Reference Validation

**Check:**
- All internal links resolve to existing files
- All prerequisite documents exist
- All related documents exist
- No circular dependencies in prerequisites
- No broken relative paths

**Process:**
- Parse all markdown files
- Extract all internal links: `[text](path.md)`, `[text](../path.md)`
- Extract prerequisite and related frontmatter fields
- Verify each reference resolves to an existing document
- Check for circular prerequisite chains

**Report:**
```
Cross-Reference Validation
==========================
Checked: 127 internal references across 30 documents

✓ All internal links valid (85/85)
✓ All prerequisites exist (42/42)
✓ All related links exist (38/38)
✓ No circular dependencies found
✓ All relative paths resolve correctly

Status: PASSED
```

**On Failure:**
```
Cross-Reference Validation
==========================
Checked: 127 internal references across 30 documents

✗ Broken internal links (3):
  • workflows/smlm-analysis.md:42 → how-to/advanced-fitting.md (not found)
  • api-reference/smi-core.md:156 → ../workflows/roi-analysis.md (not found)
  • examples/multi-channel.md:23 → troubleshooting/registration.md (not found)

✓ All prerequisites exist (42/42)
⚠ Circular dependency detected:
  • workflows/a.md → workflows/b.md → workflows/a.md

Status: FAILED
Action: Fix broken links, resolve circular dependencies
```

### 3. Frontmatter Consistency

**Check:**
- All required frontmatter fields present
- Field values are valid (e.g., category matches allowed list)
- Dates are properly formatted
- Tags are consistent
- Prerequisites and related lists are arrays
- No duplicate IDs

**Required Fields:**
- title
- category
- level
- tags (array)
- summary
- estimated_time
- last_updated (YYYY-MM-DD format)
- status

**Valid Values:**
- category: getting-started, core-concepts, workflows, how-to, api-reference, examples, troubleshooting, reference
- level: beginner, intermediate, advanced
- status: complete, in_progress, needs_review

**Report:**
```
Frontmatter Consistency
=======================
Checked: 30 documents

✓ All required fields present (30/30)
✓ All categories valid (30/30)
✓ All levels valid (30/30)
✓ All dates properly formatted (30/30)
✓ All tags are arrays (30/30)
✓ No duplicate document IDs

Status: PASSED
```

**On Failure:**
```
Frontmatter Consistency
=======================
Checked: 30 documents

✗ Missing required fields (2):
  • workflows/bagol-clustering.md: missing 'estimated_time'
  • examples/tracking.md: missing 'summary'

✗ Invalid category (1):
  • how-to/custom-processing.md: category='advanced' (not in allowed list)

✗ Invalid date format (1):
  • api-reference/smi-core.md: last_updated='15-01-2025' (should be 2025-01-15)

✓ All levels valid (30/30)
✓ All tags are arrays (30/30)

Status: FAILED
Action: Fix frontmatter in listed documents
```

### 4. Manifest Accuracy

**Check:**
- Manifest matches actual files on disk
- Document counts are accurate
- All documents listed in manifest have entries
- Category document counts match reality
- No duplicate entries
- Generated timestamp is reasonable

**Report:**
```
Manifest Accuracy
=================
✓ Document count matches (30 in manifest, 30 on disk)
✓ All manifest entries have corresponding files
✓ Category counts accurate
✓ No duplicate entries
✓ Generated timestamp valid (2025-01-15T10:30:00Z)
✓ All frontmatter synchronized with manifest

Status: PASSED
```

**On Failure:**
```
Manifest Accuracy
=================
✗ Document count mismatch: 30 in manifest, 28 on disk
✗ Missing files (2):
  • examples/3d-analysis.md (in manifest, not on disk)
  • troubleshooting/gpu.md (in manifest, not on disk)

⚠ Category count incorrect:
  • getting-started: manifest says 4, actual 3

✓ No duplicate entries
✗ Generated timestamp is 3 months old (2024-10-15T10:30:00Z)

Status: FAILED
Action: Regenerate manifest: /smite_update_llm_docs
```

### 5. Link Validation (External)

**Check:**
- External links are properly formatted
- GitHub links point to correct repository
- Documentation links are accessible (optional)

**Process:**
- Extract all external URLs
- Check URL format validity
- Verify GitHub URLs point to LidkeLab/smite
- Optionally ping URLs to check accessibility

**Report:**
```
External Link Validation
========================
Found: 45 external links across 30 documents

✓ All URLs properly formatted (45/45)
✓ All GitHub links correct (12/12)
⚠ Accessibility not tested (use --check-live to test)

Status: PASSED
```

### 6. Code Example Syntax

**Check:**
- All MATLAB code blocks have valid syntax
- Code blocks are properly formatted
- No obvious syntax errors (missing quotes, parentheses, etc.)

**Process:**
- Extract all ```matlab code blocks
- Check for:
  - Balanced parentheses, brackets, braces
  - Balanced quotes
  - Valid MATLAB keywords
  - Basic syntax structure
- Report line numbers of issues

**Report:**
```
Code Example Syntax
===================
Checked: 47 MATLAB code blocks across 30 documents

✓ All code blocks properly formatted (47/47)
✓ Balanced parentheses/brackets (47/47)
✓ Balanced quotes (47/47)
✓ No obvious syntax errors

Status: PASSED
```

**On Failure:**
```
Code Example Syntax
===================
Checked: 47 MATLAB code blocks across 30 documents

✗ Syntax errors detected (2):
  • workflows/smlm-analysis.md:125
    Error: Unbalanced parentheses
    Code: SMF = smi_core.SingleMoleculeFitting(

  • examples/basic-localization.md:67
    Error: Unterminated string
    Code: disp('Result: %f, SMD.X(1));

✓ Balanced brackets (45/47)
✓ No obvious syntax errors in other blocks (45/47)

Status: FAILED
Action: Fix syntax errors in listed code blocks
```

### 7. Code Example Execution (Optional)

**Check:**
- Code examples actually run in MATLAB
- Output matches expected output
- No runtime errors

**Note:** This is optional and computationally expensive.

**Process:**
- Extract executable code blocks (those marked for execution)
- Run in MATLAB environment
- Capture output and errors
- Compare with expected output if provided

**Report:**
```
Code Example Execution
======================
⚠ Code execution test skipped (use --execute to test)

To test code execution:
  /smite_validate_llm_docs --execute

Note: This will run all code examples and may take several minutes.

Status: SKIPPED
```

**If executed:**
```
Code Example Execution
======================
Executed: 47 code blocks across 30 documents

✓ Successful execution (45/47)
✗ Runtime errors (2):
  • workflows/bagol-clustering.md:134
    Error: Undefined variable 'ROI'

  • examples/tracking-diffusion.md:89
    Error: Input data dimensions mismatch

⚠ Output mismatch (1):
  • examples/basic-localization.md:56
    Expected: "Found 47 localizations"
    Got: "Found 48 localizations"

Status: FAILED
Action: Fix code examples with runtime errors
```

### 8. Document Structure

**Check:**
- All documents follow template structure
- Required sections present
- Consistent heading levels
- Proper markdown formatting

**Required Sections:**
- Purpose
- Prerequisites (if applicable)
- Overview
- Main content sections
- See Also
- Next Steps (if applicable)

**Report:**
```
Document Structure
==================
Checked: 30 documents

✓ All documents follow template (30/30)
✓ Required sections present (30/30)
✓ Consistent heading hierarchy (30/30)
✓ Proper markdown formatting (30/30)

Status: PASSED
```

## Validation Report Format

### Summary

```
LLM Documentation Validation Report
====================================
Validation Date: 2025-01-20 14:30:00
Documents Validated: 30
Status File: .claude/llm-docs-status.json

Overall Status: ⚠ PASSED WITH WARNINGS

Test Results:
  ✓ File Existence:            PASSED
  ✓ Cross-References:           PASSED
  ✓ Frontmatter Consistency:    PASSED
  ✓ Manifest Accuracy:          PASSED
  ✓ External Links:             PASSED
  ⚠ Code Example Syntax:        WARNING (1 minor issue)
  ○ Code Example Execution:     SKIPPED
  ✓ Document Structure:         PASSED

Health Score: 96% (Excellent)

Issues Found:
  Critical: 0
  Major:    0
  Minor:    1
  Warnings: 2

Details below...
```

### Detailed Results

Show detailed results for each check as formatted above.

### Issue Summary

```
Issue Summary
=============

Minor Issues (1):
  1. workflows/smlm-analysis.md:125
     Missing period at end of sentence

Warnings (2):
  1. Code execution not tested (15 untested code blocks)
  2. External link accessibility not verified (45 unchecked URLs)

Recommendations:
  • Fix minor formatting issue in workflows/smlm-analysis.md
  • Consider running validation with --execute flag
  • Run validation with --check-live for link checking
```

### Update Status File

After validation, update `.claude/llm-docs-status.json`:

```json
{
  "validation": {
    "last_run": "2025-01-20T14:30:00Z",
    "overall_status": "passed_with_warnings",
    "health_score": 96.0,
    "tests": {
      "file_existence": "passed",
      "cross_references": "passed",
      "frontmatter": "passed",
      "manifest": "passed",
      "external_links": "passed",
      "code_syntax": "warning",
      "code_execution": "skipped",
      "document_structure": "passed"
    },
    "issues": {
      "critical": 0,
      "major": 0,
      "minor": 1,
      "warnings": 2
    },
    "files_checked": 30,
    "cross_refs_checked": 127,
    "code_blocks_checked": 47
  }
}
```

## Command Options

Support optional flags:

- `--execute` : Run code examples in MATLAB
- `--check-live` : Verify external links are accessible
- `--fix-auto` : Automatically fix simple issues (whitespace, formatting)
- `--detailed` : Show detailed output for all checks
- `--json` : Output results in JSON format

## Exit Codes / Status

Return clear status:
- **PASSED**: All critical checks passed, no issues
- **PASSED WITH WARNINGS**: All checks passed but minor issues noted
- **FAILED**: One or more critical checks failed
- **ERROR**: Validation could not complete

## Important Notes

- Run after building or updating documentation
- Recommend running before commits
- Can be automated in CI/CD
- Some checks (code execution) are optional but valuable
- Update status file with validation results
- Provide actionable fix suggestions
- Don't fail on warnings, only on actual errors
- Be forgiving of minor formatting issues
