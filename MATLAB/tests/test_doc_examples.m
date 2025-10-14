function results = test_doc_examples()
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

    fprintf('Testing %d code examples from documentation...\n', length(blocks));
    fprintf('========================================\n\n');

    % Initialize results
    results = struct();
    results.total = length(blocks);
    results.passed = 0;
    results.failed = 0;
    results.skipped = 0;
    results.syntax_only = 0;
    results.details = cell(length(blocks), 1);

    % Test each block
    for i = 1:length(blocks)
        block = blocks(i);

        if mod(i, 100) == 0
            fprintf('[%d/%d] Progress: %d%% complete\n', i, length(blocks), ...
                    round(100*i/length(blocks)));
        end

        result = test_code_block(block, i);
        results.details{i} = result;

        switch result.status
            case 'passed'
                results.passed = results.passed + 1;
            case 'failed'
                results.failed = results.failed + 1;
                fprintf('[%d] FAILED: %s:%d\n', i, result.file, result.line);
                fprintf('  Error: %s\n', result.error);
            case 'skipped'
                results.skipped = results.skipped + 1;
            case 'syntax_only'
                results.syntax_only = results.syntax_only + 1;
        end
    end

    % Summary
    fprintf('\n========================================\n');
    fprintf('Summary\n');
    fprintf('========================================\n');
    fprintf('Total:        %d\n', results.total);
    fprintf('Passed:       %d (%.1f%%)\n', results.passed, ...
            100*results.passed/results.total);
    fprintf('Failed:       %d (%.1f%%)\n', results.failed, ...
            100*results.failed/results.total);
    fprintf('Skipped:      %d (%.1f%%)\n', results.skipped, ...
            100*results.skipped/results.total);
    fprintf('Syntax only:  %d (%.1f%%)\n', results.syntax_only, ...
            100*results.syntax_only/results.total);

    % Save results
    output_file = fullfile(fileparts(mfilename('fullpath')), ...
                           'doc_examples_results.mat');
    save(output_file, 'results');
    fprintf('\nResults saved to: %s\n', output_file);

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
    result.category = '';

    % Classify the code block
    [is_testable, category, skip_reason] = classify_code_block(block.code, block.context);

    if ~is_testable
        result.status = category;
        result.reason = skip_reason;
        result.category = category;
        return;
    end

    % Try to execute the code
    try
        % Capture warnings
        lastwarn('');
        orig_warn = warning('query', 'all');
        warning('on', 'all');

        % Execute code
        evalc_output = evalc('eval(block.code)');

        % Check for warnings
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            result.warnings = {warnId, warnMsg};
        end

        % Restore warnings
        warning(orig_warn);

        result.status = 'passed';
        result.output = evalc_output;
        result.category = 'executable';

    catch ME
        result.status = 'failed';
        result.error = sprintf('%s: %s', ME.identifier, ME.message);
        result.stack_trace = arrayfun(@(s) sprintf('%s:%d', s.name, s.line), ...
                                      ME.stack, 'UniformOutput', false);
        result.category = 'executable';
    end
end

function [is_testable, category, reason] = classify_code_block(code, context)
% Classify whether a code block is testable and what type it is
%
% Returns:
%   is_testable - boolean, whether to attempt execution
%   category - 'syntax_only', 'skipped', or 'executable'
%   reason - explanation string

    is_testable = true;
    category = 'executable';
    reason = '';

    % Remove leading/trailing whitespace
    code = strtrim(code);

    % Empty code
    if isempty(code)
        is_testable = false;
        category = 'skipped';
        reason = 'Empty code block';
        return;
    end

    % Comments only
    lines = strsplit(code, '\n');
    non_comment_lines = lines(~startsWith(strtrim(lines), '%'));
    if isempty(non_comment_lines) || all(cellfun(@isempty, strtrim(non_comment_lines)))
        is_testable = false;
        category = 'syntax_only';
        reason = 'Comments only';
        return;
    end

    % Single line syntax examples (no assignment, no function call)
    if length(non_comment_lines) == 1
        line = strtrim(non_comment_lines{1});
        % Check if it's just showing syntax: obj.property or ClassName.method
        if ~contains(line, '=') && ~contains(line, '(') && contains(line, '.')
            is_testable = false;
            category = 'syntax_only';
            reason = 'Syntax example (property access pattern)';
            return;
        end
        % Check if it's incomplete (ends with ... or just partial)
        if endsWith(line, '...') || (length(line) < 10 && ~contains(line, '()'))
            is_testable = false;
            category = 'syntax_only';
            reason = 'Incomplete statement';
            return;
        end
    end

    % Function signature examples (show input/output pattern)
    if startsWith(code, '[') && contains(code, '] = ') && length(non_comment_lines) == 1
        % e.g., [Output1, Output2] = obj.methodName(Input1, Input2)
        if contains(code, 'obj.') || ~contains(code, ';')
            is_testable = false;
            category = 'syntax_only';
            reason = 'Function signature example';
            return;
        end
    end

    % GUI-related code
    if contains(code, '.gui()') || contains(code, 'uigetfile') || ...
       contains(code, 'uiputfile') || contains(code, 'inputdlg') || ...
       contains(code, 'uifigure') || contains(code, 'uicontrol')
        is_testable = false;
        category = 'skipped';
        reason = 'Requires GUI interaction';
        return;
    end

    % Requires external files that won't exist
    if (contains(code, 'h5read') || contains(code, 'load(')) && ...
       ~contains(code, 'simulate') && ~contains(code, 'test') && ...
       ~contains(code, 'example')
        is_testable = false;
        category = 'skipped';
        reason = 'Requires external data files';
        return;
    end

    % Uses undefined variables (placeholder examples)
    % Common placeholder patterns
    placeholders = {'Input1', 'Input2', 'Output1', 'Output2', 'varargin', ...
                   'obj', 'ClassName', 'methodName', 'PropertyName', ...
                   'newValue', 'inputs', 'result', 'value'};

    % Check if code uses placeholders without defining them
    for p = 1:length(placeholders)
        placeholder = placeholders{p};
        % If placeholder appears but is never assigned
        if contains(code, placeholder) && ~contains(code, [placeholder, ' ='])
            % Exception: varargin is valid in function definitions
            if strcmp(placeholder, 'varargin') && contains(code, 'function')
                continue;
            end
            is_testable = false;
            category = 'syntax_only';
            reason = sprintf('Uses placeholder variable: %s', placeholder);
            return;
        end
    end

    % Context indicates it's documentation
    if contains(context, 'syntax') || contains(context, 'pattern') || ...
       contains(context, 'format') || contains(context, 'signature')
        % Check if code looks incomplete
        if length(non_comment_lines) <= 3 && ~any(contains(code, {'function ', 'for ', 'if '}))
            is_testable = false;
            category = 'syntax_only';
            reason = 'Syntax documentation example';
            return;
        end
    end

    % Multi-line but all are property assignments without context
    if all(contains(non_comment_lines, '.') & contains(non_comment_lines, '='))
        first_line = non_comment_lines{1};
        if contains(first_line, '.')
            parts = strsplit(first_line, '.');
            var_name = strtrim(parts{1});
            % Check if variable is defined in code
            if ~any(contains(non_comment_lines, [var_name, ' = ']))
                is_testable = false;
                category = 'skipped';
                reason = 'Property assignments without object creation';
                return;
            end
        end
    end

    % Mex/compilation commands
    if contains(code, 'mex ') || contains(code, 'nvcc') || contains(code, 'gcc')
        is_testable = false;
        category = 'skipped';
        reason = 'Compilation command';
        return;
    end

    % Installation/setup commands
    if contains(code, 'addpath') || contains(code, 'savepath') || ...
       contains(code, 'cd ') || contains(code, 'git clone')
        is_testable = false;
        category = 'skipped';
        reason = 'Installation/setup command';
        return;
    end
end

function generate_report(results)
% Generate detailed audit report

    report_file = fullfile(fileparts(mfilename('fullpath')), ...
                           'doc_examples_audit_report.txt');

    fid = fopen(report_file, 'w');

    fprintf(fid, 'LLM Documentation Audit Report\n');
    fprintf(fid, '==============================\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));

    fprintf(fid, 'Summary\n');
    fprintf(fid, '-------\n');
    fprintf(fid, 'Total code blocks: %d\n', results.total);
    fprintf(fid, 'Passed:       %d (%.1f%%)\n', results.passed, ...
            100*results.passed/results.total);
    fprintf(fid, 'Failed:       %d (%.1f%%)\n', results.failed, ...
            100*results.failed/results.total);
    fprintf(fid, 'Skipped:      %d (%.1f%%)\n', results.skipped, ...
            100*results.skipped/results.total);
    fprintf(fid, 'Syntax only:  %d (%.1f%%)\n\n', results.syntax_only, ...
            100*results.syntax_only/results.total);

    % Categorize failures
    if results.failed > 0
        fprintf(fid, '\nFailed Tests\n');
        fprintf(fid, '============\n\n');

        error_categories = struct();

        for i = 1:length(results.details)
            detail = results.details{i};
            if strcmp(detail.status, 'failed')
                % Categorize error
                err_id = detail.error;
                if contains(err_id, 'Undefined')
                    category = 'Undefined variable/function';
                elseif contains(err_id, 'Too many') || contains(err_id, 'Not enough')
                    category = 'Incorrect arguments';
                elseif contains(err_id, 'Class') || contains(err_id, 'method')
                    category = 'Class/method error';
                else
                    category = 'Other';
                end

                if ~isfield(error_categories, matlab.lang.makeValidName(category))
                    error_categories.(matlab.lang.makeValidName(category)) = {};
                end
                error_categories.(matlab.lang.makeValidName(category)){end+1} = i;
            end
        end

        % Print by category
        category_names = fieldnames(error_categories);
        for c = 1:length(category_names)
            cat_name = category_names{c};
            % Convert back to readable name
            readable_name = strrep(cat_name, '_', ' ');
            indices = error_categories.(cat_name);

            fprintf(fid, '%s (%d failures)\n', readable_name, length(indices));
            fprintf(fid, '%s\n\n', repmat('-', 1, length(readable_name)+15));

            for j = 1:min(10, length(indices))  % Show first 10
                idx = indices{j};
                detail = results.details{idx};
                fprintf(fid, '%d. %s:%d\n', idx, detail.file, detail.line);
                fprintf(fid, '   Error: %s\n', detail.error);
                fprintf(fid, '   Code:\n');
                code_lines = strsplit(detail.code, sprintf('\n'));
                for k = 1:min(5, length(code_lines))
                    fprintf(fid, '      %s\n', code_lines{k});
                end
                if length(code_lines) > 5
                    fprintf(fid, '      ... (%d more lines)\n', length(code_lines)-5);
                end
                fprintf(fid, '\n');
            end

            if length(indices) > 10
                fprintf(fid, '   ... and %d more failures in this category\n\n', ...
                        length(indices)-10);
            end
        end
    end

    % Warnings summary
    warning_count = 0;
    for i = 1:length(results.details)
        if ~isempty(results.details{i}.warnings)
            warning_count = warning_count + 1;
        end
    end

    if warning_count > 0
        fprintf(fid, '\nWarnings\n');
        fprintf(fid, '========\n');
        fprintf(fid, '%d tests generated warnings\n\n', warning_count);
    end

    % Recommendations
    fprintf(fid, '\nRecommendations\n');
    fprintf(fid, '===============\n\n');

    if results.failed > 0
        fail_rate = 100*results.failed/results.total;
        if fail_rate > 10
            fprintf(fid, 'HIGH PRIORITY: %.1f%% of code blocks failed execution\n', fail_rate);
            fprintf(fid, '  - Review failed examples for API errors\n');
            fprintf(fid, '  - Check for undefined variables in examples\n');
            fprintf(fid, '  - Verify function signatures match implementation\n\n');
        else
            fprintf(fid, 'MEDIUM PRIORITY: %.1f%% of code blocks failed\n', fail_rate);
            fprintf(fid, '  - Most failures likely due to missing test setup\n');
            fprintf(fid, '  - Review individual cases\n\n');
        end
    else
        fprintf(fid, 'EXCELLENT: No executable code blocks failed\n\n');
    end

    fprintf(fid, 'Summary:\n');
    fprintf(fid, '  - %.1f%% of blocks are executable code\n', ...
            100*(results.passed+results.failed)/results.total);
    fprintf(fid, '  - %.1f%% are syntax/documentation examples\n', ...
            100*results.syntax_only/results.total);
    fprintf(fid, '  - %.1f%% require special setup (skipped)\n', ...
            100*results.skipped/results.total);

    fclose(fid);
    fprintf('\nDetailed report saved to: %s\n', report_file);
end
