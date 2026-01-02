#!/usr/bin/env python3
"""
Parse Snakemake dry-run output and extract shell commands in execution order.
"""

import argparse
import re
import sys
from pathlib import Path


def parse_snakemake_output(input_file):
    """
    Parse snakemake dry-run output and extract shell commands.
    
    Args:
        input_file: Path to file containing snakemake -p -n --forceall output
        
    Returns:
        List of tuples (rule_name, job_number, command)
    """
    commands = []
    current_rule = None
    current_job = None
    in_shell_block = False
    shell_lines = []
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            # Match rule lines like "rule rule_name:"
            rule_match = re.match(r'^rule\s+(\S+):\s*$', line)
            if rule_match:
                # Save previous command if exists
                if in_shell_block and shell_lines and current_rule and current_job:
                    command = '\n'.join(shell_lines).strip()
                    # Replace multiple spaces with single space
                    command = re.sub(r' {2,}', ' ', command)
                    if command:
                        commands.append((current_rule, current_job, command))
                    shell_lines = []
                    in_shell_block = False
                
                current_rule = rule_match.group(1)
                current_job = None
                continue
            
            # Extract job number from output line like "    output: results/1-sample_name/"
            if current_rule:
                output_match = re.match(r'^\s+output:\s+results/(\d+)-', line)
                if output_match:
                    current_job = int(output_match.group(1))
                    continue
            
            # Detect shell command block (starts after output/input/wildcards sections)
            if current_rule and current_job is not None:
                # Skip metadata lines
                if re.match(r'^\s+(Reason|Input|Output|Wildcards|jobid|threads|resources|reason):', line):
                    continue
                
                # Shell command starts with indentation and actual command
                if re.match(r'^\s+\S', line):
                    in_shell_block = True
                    shell_lines.append(line)
                elif in_shell_block:
                    if line.strip():
                        shell_lines.append(line)
                    else:
                        # Empty line ends shell block
                        if shell_lines:
                            command = '\n'.join(shell_lines).strip()
                            # Replace multiple spaces with single space
                            command = re.sub(r' {2,}', ' ', command)
                            if command:
                                commands.append((current_rule, current_job, command))
                            shell_lines = []
                        in_shell_block = False
    
    # Handle last command
    if in_shell_block and shell_lines and current_rule and current_job:
        command = '\n'.join(shell_lines).strip()
        # Replace multiple spaces with single space
        command = re.sub(r' {2,}', ' ', command)
        if command:
            commands.append((current_rule, current_job, command))
    
    return commands


def write_output(commands, output_file):
    """
    Write extracted commands to output file.
    
    Args:
        commands: List of tuples (rule_name, job_number, command)
        output_file: Path to output file
    """
    with open(output_file, 'w') as f:
        for rule_name, job_num, command in commands:
            f.write(f"# Job {job_num}: {rule_name}\n")
            f.write(f"{command}\n")
            f.write("-" * 80 + "\n\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract shell commands from Snakemake dry-run output in execution order"
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file containing snakemake -p -n --forceall output'
    )
    parser.add_argument(
        '-o', '--output',
        default='snakemake_commands.txt',
        help='Output file (default: snakemake_commands.txt)'
    )
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = Path(args.output).parent
    if output_dir != Path('.'):
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse and write output
    print(f"Parsing {args.input}...")
    commands = parse_snakemake_output(args.input)
    print(f"Found {len(commands)} shell commands")
    
    write_output(commands, args.output)
    print(f"Commands written to {args.output}")


if __name__ == "__main__":
    main()