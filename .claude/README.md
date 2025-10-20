# .claude Directory

This directory contains Claude Code-specific configuration and tracking files for the smite project.

## Contents

### Slash Commands (`commands/`)

Custom slash commands for managing LLM documentation:

- **`/smite_build_llm_docs`** - Build LLM documentation in phases
  - Creates modular, LLM-friendly documentation in `doc/llm-guide/`
  - Builds in 3 phases: Essential → Comprehensive → Maintenance
  - Tracks progress and enables resumable builds
  - See: `commands/smite_build_llm_docs.md`

- **`/smite_update_llm_docs`** - Update documentation when code changes
  - Analyzes git changes since last update
  - Maps changed files to affected documentation
  - Proposes and executes targeted updates
  - Maintains consistency across documentation
  - See: `commands/smite_update_llm_docs.md`

- **`/smite_docs_status`** - Quick status report
  - Shows build progress by phase
  - Lists pending updates
  - Reports documentation health
  - Provides actionable recommendations
  - See: `commands/smite_docs_status.md`

- **`/smite_validate_llm_docs`** - Validate documentation quality
  - Checks file existence and structure
  - Validates cross-references and links
  - Verifies frontmatter consistency
  - Tests code example syntax
  - Optionally executes code examples
  - See: `commands/smite_validate_llm_docs.md`

### Status Tracking

- **`llm-docs-status.json`** - Documentation build state (created on first use)
  - Tracks which documents are complete
  - Records pending updates
  - Stores validation results
  - Enables resumable builds
  - Schema: `llm-docs-status.schema.json`

- **`llm-docs-status.schema.json`** - JSON schema for status file
  - Defines structure and validation rules
  - Documents all fields and their meanings
  - Reference for tooling and automation

## Quick Start

### Building Documentation

First time:
```
/smite_build_llm_docs
```

This will:
1. Create `doc/llm-guide/` directory structure
2. Build Phase 1 (12 essential documents)
3. Generate `manifest.json` and `index.md`
4. Create `llm-docs-status.json` to track progress

Continue building:
```
/smite_build_llm_docs
```

Will automatically continue to Phase 2.

### Maintaining Documentation

After making code changes:
```
/smite_update_llm_docs
```

This will:
1. Analyze git changes since last update
2. Identify affected documentation
3. Propose and apply updates
4. Keep documentation in sync with code

### Checking Status

Anytime:
```
/smite_docs_status
```

Shows:
- Build progress
- Pending updates
- Documentation health
- Next recommended actions

### Validating Quality

Before commits or releases:
```
/smite_validate_llm_docs
```

Validates:
- All files exist
- Cross-references work
- Frontmatter is consistent
- Code examples are syntactically correct
- Links are valid

## Workflow Integration

### Recommended Workflow

1. **Initial Setup** (one time):
   ```
   /smite_build_llm_docs
   ```

2. **During Development**:
   - Make code changes
   - Commit changes
   - Update docs: `/smite_update_llm_docs`
   - Validate: `/smite_validate_llm_docs`

3. **Periodic Maintenance**:
   - Check status: `/smite_docs_status`
   - Build more docs: `/smite_build_llm_docs`
   - Review pending updates

4. **Before Releases**:
   - Update all docs: `/smite_update_llm_docs`
   - Validate: `/smite_validate_llm_docs`
   - Review coverage: `/smite_docs_status`

### Git Integration

The `llm-docs-status.json` file can be:

**Committed (Recommended):**
- ✅ Team visibility into doc progress
- ✅ Track documentation state over time
- ✅ Enable collaboration
- ❌ May cause merge conflicts (rare)

**Gitignored:**
- ✅ No merge conflicts
- ✅ Personal tracking only
- ❌ No team visibility
- ❌ Lose history

We recommend **committing** it for team visibility. Treat it like `package-lock.json`.

### CI/CD Integration (Future)

These commands could be automated:

```yaml
# .github/workflows/llm-docs.yml
name: LLM Docs Check

on:
  pull_request:
    paths:
      - 'MATLAB/**'
      - 'doc/**'

jobs:
  check-docs:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Check if docs need updates
        run: |
          # Parse llm-docs-status.json
          # Check if PR changes require doc updates
          # Comment on PR if updates needed
```

## Documentation Structure

The generated documentation lives in `doc/llm-guide/`:

```
doc/llm-guide/
  index.md                   # Master navigation
  manifest.json              # Machine-readable catalog

  getting-started/           # Beginner guides
  core-concepts/             # Foundational knowledge
  workflows/                 # Complete analysis pipelines
  how-to/                    # Task-based guides
  api-reference/             # Namespace documentation
  examples/                  # Practical examples
  troubleshooting/           # Problem-solving
  reference/                 # Technical details
```

Each document includes:
- Complete frontmatter metadata
- Working code examples
- Cross-references to related docs
- Prerequisites and next steps
- Troubleshooting tips

## Access Methods

Documentation can be accessed by LLMs via:

1. **Local file paths** (during development)
   ```
   file:///C:/Users/klidke/Documents/MATLAB/smite/doc/llm-guide/
   ```

2. **GitHub Raw URLs** (after pushing)
   ```
   https://raw.githubusercontent.com/LidkeLab/smite/main/doc/llm-guide/
   ```

3. **GitHub Pages** (if enabled)
   ```
   https://lidkelab.github.io/smite/llm-guide/
   ```

4. **MCP Server** (future - custom implementation)
   - Provides semantic search
   - Dynamic document retrieval
   - Context-aware suggestions

## Status File Structure

Example `llm-docs-status.json`:

```json
{
  "version": "1.0.0",
  "base_path": "doc/llm-guide/",
  "last_build": "2025-01-10T14:30:00Z",
  "last_update": "2025-01-15T09:00:00Z",

  "phases": {
    "phase1_essential": {
      "status": "complete",
      "completed": "2025-01-10T14:30:00Z",
      "target_docs": 12,
      "completed_docs": 12
    },
    "phase2_comprehensive": {
      "status": "in_progress",
      "started": "2025-01-12T09:00:00Z",
      "target_docs": 40,
      "completed_docs": 18
    }
  },

  "documents": {
    "getting-started/quickstart.md": {
      "status": "complete",
      "created": "2025-01-10T14:00:00Z",
      "updated": "2025-01-10T14:00:00Z",
      "word_count": 650,
      "validation": "passed",
      "source_files": ["README.md", "doc/CoreOverview.md"]
    }
  },

  "pending_updates": [],

  "coverage": {
    "total_planned": 52,
    "total_complete": 30,
    "percentage": 57.7
  }
}
```

## Design Philosophy

### Why Slash Commands?

- **Self-Documenting**: Process is encoded in the commands themselves
- **Discoverable**: Easy to find with `ls .claude/commands`
- **Maintainable**: Update commands as needs evolve
- **Repeatable**: Same process every time
- **Collaborative**: Team uses same workflow

### Why Status Tracking?

- **Resumable**: Can pause and continue builds
- **Transparent**: See exactly what's done and what's pending
- **Intelligent**: Updates only what changed
- **Validated**: Track documentation health
- **Historical**: See progress over time

### Why Modular Documentation?

- **LLM-Friendly**: Load only relevant pieces
- **Maintainable**: Update individual docs independently
- **Discoverable**: Manifest enables semantic search
- **Progressive**: Beginner → Advanced learning paths
- **Task-Oriented**: "How do I..." gets quick answers

## Troubleshooting

### Status file corrupted
Delete `.claude/llm-docs-status.json` and rebuild:
```
rm .claude/llm-docs-status.json
/smite_build_llm_docs
```

### Commands not found
Make sure you're in the smite repository root:
```
cd ~/Documents/MATLAB/smite
/smite_docs_status
```

### Documentation out of sync
Force update all docs:
```
/smite_update_llm_docs
```

### Validation failures
Fix issues and re-validate:
```
/smite_validate_llm_docs
```

## Future Enhancements

Potential additions:
- [ ] Automated documentation generation from docstrings
- [ ] MCP server for semantic search
- [ ] CI/CD integration for automatic updates
- [ ] Version-specific documentation branches
- [ ] Interactive documentation testing
- [ ] Documentation metrics and analytics
- [ ] Automatic issue detection and suggestions

## Contributing

To add new commands:
1. Create `.claude/commands/new_command.md`
2. Follow existing command structure
3. Document in this README
4. Test thoroughly
5. Update if needed

## Resources

- [Claude Code Documentation](https://docs.claude.com/claude-code)
- [Slash Commands Guide](https://docs.claude.com/claude-code/guides/slash-commands)
- [smite Main README](../README.md)
- [LLM Guide Index](../doc/llm-guide/index.md) (after building)
