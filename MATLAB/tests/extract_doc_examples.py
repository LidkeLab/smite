#!/usr/bin/env python3
"""
Extract MATLAB code examples from LLM documentation for testing.

This script scans all markdown files in doc/llm-guide/ and extracts
code blocks marked with ```matlab. Each code block is saved with
metadata about its source file and line number.
"""

import re
import json
from pathlib import Path
from typing import List, Dict, Tuple

def extract_code_blocks(md_file: Path) -> List[Dict]:
    """Extract all ```matlab code blocks from a markdown file.

    Returns list of dicts with:
    - file: source file path (relative to repo root)
    - line: starting line number
    - code: the code block content
    - context: surrounding text for context
    """
    with open(md_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    blocks = []
    in_code_block = False
    current_block = []
    start_line = 0
    block_lang = None

    for i, line in enumerate(lines, 1):
        # Check for code block start
        match = re.match(r'^```(\w+)?', line)
        if match and not in_code_block:
            block_lang = match.group(1) or 'unknown'
            if block_lang.lower() == 'matlab':
                in_code_block = True
                current_block = []
                start_line = i + 1
        elif line.startswith('```') and in_code_block:
            # End of code block
            code = ''.join(current_block)

            # Get context (lines before code block)
            context_start = max(0, start_line - 10)
            context = ''.join(lines[context_start:start_line-2])

            blocks.append({
                'file': str(md_file.relative_to(Path.cwd())),
                'line': start_line,
                'code': code,
                'context': context.strip()
            })

            in_code_block = False
            current_block = []
        elif in_code_block:
            current_block.append(line)

    return blocks

def main():
    """Extract all code blocks from documentation."""
    repo_root = Path.cwd()
    doc_dir = repo_root / 'doc' / 'llm-guide'

    if not doc_dir.exists():
        print(f"Error: Documentation directory not found: {doc_dir}")
        return

    all_blocks = []
    file_count = 0

    # Find all markdown files
    md_files = sorted(doc_dir.rglob('*.md'))

    print(f"Scanning {len(md_files)} markdown files...")

    for md_file in md_files:
        blocks = extract_code_blocks(md_file)
        if blocks:
            file_count += 1
            all_blocks.extend(blocks)
            print(f"  {md_file.relative_to(doc_dir)}: {len(blocks)} code blocks")

    print(f"\nTotal: {len(all_blocks)} code blocks from {file_count} files")

    # Save to JSON
    output_file = repo_root / 'MATLAB' / 'tests' / 'doc_examples.json'
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(all_blocks, f, indent=2)

    print(f"\nSaved to: {output_file}")

    # Generate MATLAB test script
    generate_matlab_test_script(all_blocks, repo_root)

def generate_matlab_test_script(blocks: List[Dict], repo_root: Path):
    """Generate MATLAB script to test all code blocks."""

    output_file = repo_root / 'MATLAB' / 'tests' / 'test_doc_examples.m'

    script = """function results = test_doc_examples()
% Test all MATLAB code examples from LLM documentation
%
% This script extracts code blocks from doc_examples.json and
% executes each one in an isolated environment, capturing errors
% and warnings.
%
% Results are saved to doc_examples_results.mat

    % Load extracted code blocks
    json_file = fullfile(fileparts(mfilename('fullpath')), 'doc_examples.json');
    json_text = fileread(json_file);
    blocks = jsondecode(json_text);

    fprintf('Testing %d code examples from documentation...\\n', length(blocks));
    fprintf('========================================\\n\\n');

    % Initialize results
    results = struct();
    results.total = length(blocks);
    results.passed = 0;
    results.failed = 0;
    results.skipped = 0;
    results.details = cell(length(blocks), 1);

    % Test each block
    for i = 1:length(blocks)
        block = blocks(i);

        fprintf('[%d/%d] Testing: %s:%d\\n', i, length(blocks), ...
                block.file, block.line);

        result = test_code_block(block, i);
        results.details{i} = result;

        if strcmp(result.status, 'passed')
            results.passed = results.passed + 1;
            fprintf('  ✓ PASSED\\n\\n');
        elseif strcmp(result.status, 'failed')
            results.failed = results.failed + 1;
            fprintf('  ✗ FAILED: %s\\n\\n', result.error);
        elseif strcmp(result.status, 'skipped')
            results.skipped = results.skipped + 1;
            fprintf('  ⊘ SKIPPED: %s\\n\\n', result.reason);
        end
    end

    % Summary
    fprintf('\\n========================================\\n');
    fprintf('Summary\\n');
    fprintf('========================================\\n');
    fprintf('Total:   %d\\n', results.total);
    fprintf('Passed:  %d (%.1f%%)\\n', results.passed, ...
            100*results.passed/results.total);
    fprintf('Failed:  %d (%.1f%%)\\n', results.failed, ...
            100*results.failed/results.total);
    fprintf('Skipped: %d (%.1f%%)\\n', results.skipped, ...
            100*results.skipped/results.total);

    % Save results
    output_file = fullfile(fileparts(mfilename('fullpath')), ...
                           'doc_examples_results.mat');
    save(output_file, 'results');
    fprintf('\\nResults saved to: %s\\n', output_file);

    % Generate report
    generate_report(results);
end

function result = test_code_block(block, block_num)
% Test a single code block

    result = struct();
    result.file = block.file;
    result.line = block.line;
    result.code = block.code;
    result.status = 'unknown';
    result.error = '';
    result.warnings = {};
    result.output = '';

    % Check if code should be skipped
    skip_reason = should_skip(block.code, block.context);
    if ~isempty(skip_reason)
        result.status = 'skipped';
        result.reason = skip_reason;
        return;
    end

    % Create isolated workspace
    test_dir = fullfile(tempdir, 'smite_doc_test', sprintf('block_%04d', block_num));
    if ~exist(test_dir, 'dir')
        mkdir(test_dir);
    end
    original_dir = pwd;

    try
        % Move to test directory
        cd(test_dir);

        % Capture warnings
        warning_state = warning('query', 'all');
        warning('on', 'all');
        lastwarn('');

        % Execute code in isolated function
        output = evalc(['test_wrapper_', num2str(block_num), '(block.code)']);

        % Check for warnings
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            result.warnings = {warnId, warnMsg};
        end

        % Restore warnings
        warning(warning_state);

        result.status = 'passed';
        result.output = output;

    catch ME
        result.status = 'failed';
        result.error = sprintf('%s: %s', ME.identifier, ME.message);
        result.stack = {ME.stack.name};
    end

    % Clean up
    cd(original_dir);
    try
        rmdir(test_dir, 's');
    catch
        % Ignore cleanup errors
    end
end

function skip_reason = should_skip(code, context)
% Determine if code block should be skipped

    skip_reason = '';

    % Skip if it's just property/field examples
    if contains(code, '...') && ~contains(code, ';')
        skip_reason = 'Documentation example (not executable)';
        return;
    end

    % Skip if it requires GUI interaction
    if contains(code, '.gui()') || contains(code, 'uigetfile') || ...
       contains(code, 'uiputfile') || contains(code, 'inputdlg')
        skip_reason = 'Requires GUI interaction';
        return;
    end

    % Skip if it requires specific data files that don't exist
    if contains(code, 'load(') || contains(code, 'h5read')
        % Check if it's using example/test data
        if ~contains(code, 'test') && ~contains(code, 'example') && ...
           ~contains(code, 'simulate')
            skip_reason = 'Requires external data files';
            return;
        end
    end

    % Skip if context indicates it's a fragment
    if contains(context, 'fragment') || contains(context, 'excerpt') || ...
       contains(context, 'portion')
        skip_reason = 'Code fragment (incomplete)';
        return;
    end

    % Skip if it's just showing syntax
    lines = strsplit(strtrim(code), '\\n');
    if length(lines) <= 2 && ~contains(code, '=')
        skip_reason = 'Syntax example only';
        return;
    end
end

function generate_report(results)
% Generate detailed audit report

    report_file = fullfile(fileparts(mfilename('fullpath')), ...
                           'doc_examples_audit_report.txt');

    fid = fopen(report_file, 'w');

    fprintf(fid, 'LLM Documentation Audit Report\\n');
    fprintf(fid, '==============================\\n');
    fprintf(fid, 'Generated: %s\\n\\n', datestr(now));

    fprintf(fid, 'Summary\\n');
    fprintf(fid, '-------\\n');
    fprintf(fid, 'Total code blocks: %d\\n', results.total);
    fprintf(fid, 'Passed:  %d (%.1f%%)\\n', results.passed, ...
            100*results.passed/results.total);
    fprintf(fid, 'Failed:  %d (%.1f%%)\\n', results.failed, ...
            100*results.failed/results.total);
    fprintf(fid, 'Skipped: %d (%.1f%%)\\n\\n', results.skipped, ...
            100*results.skipped/results.total);

    % Failed tests
    if results.failed > 0
        fprintf(fid, 'Failed Tests\\n');
        fprintf(fid, '------------\\n\\n');

        for i = 1:length(results.details)
            detail = results.details{i};
            if strcmp(detail.status, 'failed')
                fprintf(fid, '%d. %s:%d\\n', i, detail.file, detail.line);
                fprintf(fid, '   Error: %s\\n', detail.error);
                fprintf(fid, '   Code:\\n');
                code_lines = strsplit(detail.code, '\\n');
                for j = 1:length(code_lines)
                    fprintf(fid, '      %s\\n', code_lines{j});
                end
                fprintf(fid, '\\n');
            end
        end
    end

    % Warnings
    warning_count = 0;
    for i = 1:length(results.details)
        if ~isempty(results.details{i}.warnings)
            warning_count = warning_count + 1;
        end
    end

    if warning_count > 0
        fprintf(fid, '\\nWarnings Generated\\n');
        fprintf(fid, '------------------\\n');
        fprintf(fid, '%d tests generated warnings\\n\\n', warning_count);
    end

    fclose(fid);
    fprintf('\\nDetailed report saved to: %s\\n', report_file);
end

% Generate wrapper functions for each block
function test_wrapper_1(code)
    eval(code);
end
"""

    # Add wrapper functions for each block (MATLAB doesn't allow dynamic function creation)
    # We'll use a single wrapper that evaluates the code

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(script)

    print(f"Generated MATLAB test script: {output_file}")

if __name__ == '__main__':
    main()
