# type: ignore

import re
import os
import sys

def parse_doxygen_comment(comment_text):
    """Parse Doxygen tags and return structured data."""
    lines = comment_text.split('\n')

    brief = ""
    description = ""
    params = []
    returns = ""
    notes = []
    warnings = []

    current_section = "description"
    current_content = []

    for line in lines:
        line = line.strip()

        # Skip empty lines and asterisks
        if not line or line == '*':
            if current_section == "description" and current_content:
                description += ' '.join(current_content) + '\n'
                current_content = []
            continue

        # Remove leading * if present
        if line.startswith('*'):
            line = line[1:].strip()

        # Parse tags
        if line.startswith('@brief'):
            brief = line.replace('@brief', '').strip()
            current_section = "brief"
        elif line.startswith('@param'):
            match = re.match(r'@param\s+(\w+)\s+(.*)', line)
            if match:
                params.append({
                    'name': match.group(1),
                    'desc': match.group(2)
                })
            current_section = "param"
        elif line.startswith('@return'):
            returns = line.replace('@return', '').strip()
            current_section = "return"
        elif line.startswith('@note'):
            notes.append(line.replace('@note', '').strip())
            current_section = "note"
        elif line.startswith('@warning'):
            warnings.append(line.replace('@warning', '').strip())
            current_section = "warning"
        elif line.startswith('@'):
            # Handle other tags like @see, @deprecated, etc.
            current_section = "other"
        else:
            # Append to current section
            if current_section == "description" and not brief:
                current_content.append(line)
            elif current_section == "return":
                returns += ' ' + line
            elif current_section == "note":
                if notes:
                    notes[-1] += ' ' + line
            elif current_section == "warning":
                if warnings:
                    warnings[-1] += ' ' + line

    # Flush remaining content
    if current_content:
        description = ' '.join(current_content)

    return {
        'brief': brief,
        'description': description.strip(),
        'params': params,
        'returns': returns.strip(),
        'notes': notes,
        'warnings': warnings
    }

def extract_public_function_signatures(header_file):
    """Extract function signatures declared in header file."""
    with open(header_file, 'r') as f:
        content = f.read()

    # Match function declarations (ending with ;)
    # Excludes typedef, struct, etc.
    pattern = r'^([a-zA-Z_][a-zA-Z0-9_*\s]*\([^)]*\))[^;{]*;'
    matches = re.findall(pattern, content, re.MULTILINE)

    signatures = set()
    for match in matches:
        sig = match.strip()
        # Filter out non-function declarations
        if '(' in sig and ')' in sig:
            signatures.add(sig)

    return signatures

def extract_function_with_signature(c_file, target_signatures):
    """Extract Doxygen comments for functions that match target signatures."""
    with open(c_file, 'r') as f:
        content = f.read()

    # Match /** ... */ followed by a function signature
    pattern = r'/\*\*(.*?)\*/(.*?)^([a-zA-Z_][a-zA-Z0-9_*\s]*\([^)]*\))[^;]*{'
    matches = re.findall(pattern, content, re.DOTALL | re.MULTILINE)

    functions = []
    for comment, spacing, signature in matches:
        sig = signature.strip()

        # Check if this signature is in the header file
        if any(sig in target_sig or target_sig in sig for target_sig in target_signatures):
            parsed = parse_doxygen_comment(comment)
            functions.append({
                'signature': sig,
                'brief': parsed['brief'],
                'description': parsed['description'],
                'params': parsed['params'],
                'returns': parsed['returns'],
                'notes': parsed['notes'],
                'warnings': parsed['warnings']
            })

    return functions

def format_function_doc(func):
    """Format a function's documentation as Markdown."""
    md = f"### `{func['signature']}`\n\n"

    if func['brief']:
        md += f"**{func['brief']}**\n\n"

    if func['description']:
        md += f"{func['description']}\n\n"

    if func['params']:
        md += "**Parameters:**\n"
        for param in func['params']:
            md += f"- `{param['name']}` — {param['desc']}\n"
        md += "\n"

    if func['returns']:
        md += f"**Returns:** {func['returns']}\n\n"

    if func['warnings']:
        for warning in func['warnings']:
            md += f"⚠️ **Warning:** {warning}\n"
        md += "\n"

    if func['notes']:
        for note in func['notes']:
            md += f"📝 **Note:** {note}\n"
        md += "\n"

    return md

def generate_readme(source_dir, output_file='API.md'):
    """Generate README from public functions (declared in headers)."""
    with open(output_file, 'w') as md:
        md.write("# API Documentation\n\n")

        for filename in sorted(os.listdir(source_dir)):
            if filename.endswith('.h'):
                header_path = os.path.join(source_dir, filename)
                c_filename = filename.replace('.h', '.c')
                c_path = os.path.join(source_dir, c_filename)

                # Get public function signatures from header
                public_sigs = extract_public_function_signatures(header_path)

                if not public_sigs or not os.path.exists(c_path):
                    continue

                # Extract documented functions from corresponding .c file
                functions = extract_function_with_signature(c_path, public_sigs)

                if functions:
                    md.write(f"## {filename}\n\n")

                    for func in functions:
                        md.write(format_function_doc(func))

if __name__ == '__main__':
    source_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    generate_readme(source_dir)
