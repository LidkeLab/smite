# Update LLM Documentation

You are tasked with updating existing LLM documentation for the smite MATLAB package when source code or documentation changes.

## Overview

This command performs incremental updates to documentation in `doc/llm-guide/` by:
- Analyzing git changes since last update
- Mapping changed files to affected documentation
- Proposing and executing targeted updates
- Maintaining documentation consistency

## Check Prerequisites

1. Verify `.claude/llm-docs-status.json` exists
   - If missing, suggest running `/smite_build_llm_docs` first
   - Cannot update what doesn't exist
2. Read status file to understand current state
3. Report current documentation status to user

## Update Detection Strategy

### Analyze Git Changes

Check what's changed since last update:

```bash
# Get last update timestamp from status file
# Then check git log since that time
git log --since="YYYY-MM-DD HH:MM:SS" --name-only --oneline --no-merges

# Also check current working tree status
git status --short
```

### Map Changes to Documentation

Use the `source_files` field in status file to identify affected docs:

**Mapping Rules:**
- `MATLAB/+smi/@SMLM/*.m` → affects:
  - `workflows/smlm-analysis.md`
  - `api-reference/smi-namespace.md`
  - `examples/basic-localization.md`

- `MATLAB/+smi/@SPT/*.m` → affects:
  - `workflows/spt-tracking.md`
  - `api-reference/smi-namespace.md`
  - `examples/tracking-diffusion.md`

- `MATLAB/+smi/@BaGoL/*.m` → affects:
  - `workflows/bagol-clustering.md`
  - `api-reference/smi-namespace.md`

- `MATLAB/+smi_core/@LocalizeData/*.m` → affects:
  - `how-to/localize-molecules.md`
  - `api-reference/smi-core.md`
  - `workflows/smlm-analysis.md`

- `doc/DataStructures/SMF.md` → affects:
  - `core-concepts/smf-structure.md`
  - `workflows/*.md` (all workflows)

- `doc/DataStructures/SMD.md` → affects:
  - `core-concepts/smd-structure.md`
  - `workflows/*.md` (all workflows)

- `MATLAB/examples/*.m` → affects:
  - Corresponding files in `examples/`

- `README.md` → affects:
  - `getting-started/quickstart.md`
  - `getting-started/installation.md`
  - `index.md`

### Classify Update Severity

For each changed file, determine impact:

**CRITICAL (Must update immediately):**
- Public API signature changes (function parameters)
- Workflow step changes
- Breaking changes
- New required parameters

**MAJOR (Should update soon):**
- New public methods/classes
- New features
- Parameter default changes
- Significant behavior changes

**MINOR (Can defer):**
- Internal refactoring
- Comment improvements
- Performance improvements
- Bug fixes that don't change usage

**IGNORE:**
- Whitespace changes
- Internal private method changes
- Test file changes (unless they reveal API changes)

## Execution Steps

### 1. Analyze Changes

Present analysis to user:
```
Documentation Update Analysis
=============================

Last Update: 2025-01-15 09:00
Current Time: 2025-01-20 14:30

Git Changes Since Last Update:
  5 commits, 8 files changed

Changed Files:
  • MATLAB/+smi/@BaGoL/combineBaGoLROIs.m (3 commits)
  • MATLAB/+smi/@BaGoL/analyzeROI.m (1 commit)
  • MATLAB/+smi/@SMLM/SMLM.m (1 commit)
  • doc/DataStructures/SMF.md (2 commits)

Affected Documentation:
  CRITICAL - Must update:
    1. workflows/bagol-clustering.md
       Reason: Core method combineBaGoLROIs changed signature
       Severity: CRITICAL (API change)
       Commits: 0a94de5, e658a5a, 999a9d2

  MAJOR - Should update:
    2. api-reference/smi-namespace.md
       Reason: BaGoL methods changed
       Severity: MAJOR (new functionality)
       Commits: 0a94de5

    3. core-concepts/smf-structure.md
       Reason: SMF structure documentation updated
       Severity: MAJOR (parameter additions)
       Commits: [commit hash]

  MINOR - Can defer:
    4. workflows/smlm-analysis.md
       Reason: Internal SMLM refactoring
       Severity: MINOR (no API impact)
       Commits: [commit hash]

Previous Pending Updates:
  • workflows/spt-tracking.md (detected 2025-01-18)

Total Documents Needing Updates: 5
```

### 2. Get User Approval

Ask user which updates to apply:
- Offer to apply all CRITICAL updates
- Suggest applying MAJOR updates
- Mention MINOR updates can be deferred
- Allow user to select specific documents
- Allow user to skip and add to pending

### 3. Execute Updates

For each document being updated:

a. **Read Current Document**
   - Parse frontmatter
   - Extract existing content
   - Note custom additions not from sources

b. **Read Updated Source Materials**
   - Read all source files listed in status
   - Read git diffs to understand exact changes
   - Check for new methods, parameters, examples

c. **Determine Update Strategy**
   - **Minor update:** Fix specific sections, preserve rest
   - **Major update:** Regenerate affected sections
   - **Full rewrite:** Regenerate entire document (rare)

d. **Generate Updated Content**
   - Follow same template as build command
   - Preserve document structure
   - Update changed sections
   - Add new sections if needed
   - Maintain code example quality

e. **Update Frontmatter**
   - Increment version or update last_updated date
   - Add to related/prerequisites if new docs exist
   - Update tags if content changed significantly

f. **Save Document**

g. **Update Status File**
   - Update `last_updated` timestamp
   - Update `word_count`
   - Add new commits to `derived_from_commits`
   - Set `validation_status` to "needs_review"
   - Update `update_reason` to "updated"

h. **Report Progress**
   ```
   [1/5] workflows/bagol-clustering.md... ✓ (1,450 words, +150 from previous)
   ```

### 4. Update Metadata Files

**If document structure changed:**
- Regenerate `manifest.json`
- Update `index.md` if new categories/docs added
- Update cross-references in affected documents

**Always:**
- Update manifest's `generated` timestamp
- Update status file's `last_update` timestamp
- Clear completed items from `pending_updates`

### 5. Detect New Pending Updates

While processing:
- Note any additional changes detected
- Add to `pending_updates` in status file
- Report to user

### 6. Final Report

```
Documentation Update Complete
=============================

Updated Documents: 4
  ✓ workflows/bagol-clustering.md (CRITICAL update)
  ✓ api-reference/smi-namespace.md (MAJOR update)
  ✓ core-concepts/smf-structure.md (MAJOR update)
  ✓ workflows/spt-tracking.md (pending update from 2025-01-18)

Deferred Updates: 1
  • workflows/smlm-analysis.md (MINOR - can update later)

New Pending Updates Detected: 0

Metadata Updated:
  ✓ manifest.json regenerated
  ✓ Status file updated

Documentation Status:
  All critical updates applied
  Documentation is current as of: 2025-01-20 14:30

Next Steps:
  • Run /smite_validate_llm_docs to verify updates
  • Review updated documents for accuracy
  • Consider applying deferred MINOR updates: /smite_update_llm_docs
  • Check status anytime: /smite_docs_status
```

## Update Strategies by Document Type

### Workflows (workflows/*.md)
- Focus on changed steps in the workflow
- Update code examples if API changed
- Update parameter explanations
- Preserve narrative flow

### How-To Guides (how-to/*.md)
- Update specific tasks affected
- Regenerate code examples
- Update parameter lists
- Check prerequisites still valid

### API Reference (api-reference/*.md)
- Add new classes/methods
- Update signatures
- Update parameter descriptions
- Add new examples

### Core Concepts (core-concepts/*.md)
- Update conceptual explanations
- Update structure diagrams if changed
- Update code examples
- Maintain pedagogical flow

### Examples (examples/*.md)
- Regenerate if source example changed
- Update if API changed
- Ensure code still runs
- Update expected output

## Smart Update Logic

### Detect API Changes
Look for in git diffs:
- Function signature changes
- New/removed parameters
- New/removed methods
- Property changes in SMF/SMD structures

### Detect Workflow Changes
Look for:
- New analysis steps
- Changed method calls
- New prerequisites
- Different output formats

### Preserve Custom Content
Do not overwrite:
- User-added notes or tips
- Custom examples not from source
- Editorial improvements
- Clarifications

Mark these sections with HTML comments if needed:
```html
<!-- CUSTOM CONTENT - DO NOT AUTO-UPDATE -->
...content...
<!-- END CUSTOM CONTENT -->
```

## Handling Edge Cases

### No Changes Detected
```
No documentation updates needed.
All docs are current as of: 2025-01-20 14:30
Last update was: 2025-01-15 09:00

Tip: Documentation is up to date!
```

### Massive Changes
If >10 documents need updates:
```
Large number of updates detected (15 documents).

Recommendations:
1. Review the changes first: git log --since="..."
2. Update in batches:
   - Update CRITICAL docs now (3 docs)
   - Update MAJOR docs next session (7 docs)
   - Defer MINOR docs (5 docs)
3. Consider if this warrants a new documentation version

Proceed with CRITICAL updates only? [Y/n]
```

### Conflicting Changes
If document was manually edited since last build:
```
Warning: workflows/smlm-analysis.md was manually edited.
Last automated update: 2025-01-15 09:00
File modified: 2025-01-18 16:30

Options:
1. Merge changes (review and combine)
2. Overwrite with generated version (lose manual edits)
3. Skip this document
4. Show diff first

Choose option [1]:
```

## Integration with Development Workflow

### Suggested Usage

**After feature work:**
```bash
git commit -m "Add new BaGoL clustering method"
# Then:
/smite_update_llm_docs
```

**Before releases:**
```bash
# Ensure all docs are current
/smite_update_llm_docs
/smite_validate_llm_docs
```

**Weekly maintenance:**
```bash
# Check for pending updates
/smite_docs_status
# Apply if needed
/smite_update_llm_docs
```

## Important Notes

- Always read status file first
- Use git to understand what changed
- Classify updates by severity
- Get user approval before major changes
- Update status file after each document
- Regenerate manifest if structure changed
- Report clearly what was done
- Add new pending items if detected
- Can be run multiple times safely
- Preserves user customizations

## Error Handling

If errors occur:
- Report what failed and why
- Save successful updates to status
- Don't roll back completed updates
- Suggest how to fix and resume
- Mark problematic docs in status file
