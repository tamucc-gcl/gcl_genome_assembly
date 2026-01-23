process FCS_CLEAN_GENOME {
  tag "${haplotype_id}"
  label 'fcs' 
  
  publishDir "${params.outdir}/contig/decontam", mode: params.publish_dir_mode

  input:
    tuple val(haplotype_id), path(assembly_fa), path(action_report)

  output:
    tuple val(haplotype_id), path("${haplotype_id}.decontaminated.fasta"), emit: decontaminated_fasta
    tuple val(haplotype_id), path("${haplotype_id}.contaminants.fasta"),   emit: contaminants_fasta
    path "${haplotype_id}.clean_stdout.log", emit: stdout_log

  script:
  """
  set -euo pipefail

  # List available scripts for debugging
  echo "Available FCS scripts:" > ${haplotype_id}.clean_stdout.log
  ls -la /app/bin/ >> ${haplotype_id}.clean_stdout.log 2>&1 || true
  echo "---" >> ${haplotype_id}.clean_stdout.log

  # Try the official script with various possible names
  CLEANED=false
  
  for script in /app/bin/clean_fasta.py /app/bin/fcs_clean.py /app/bin/apply_fcs_gx.py; do
    if [ -f "\${script}" ]; then
      echo "Using script: \${script}" | tee -a ${haplotype_id}.clean_stdout.log
      python3 "\${script}" \\
        --action-report ${action_report} \\
        --input ${assembly_fa} \\
        --output ${haplotype_id}.decontaminated.fasta \\
        --contam-output ${haplotype_id}.contaminants.fasta \\
        2>&1 | tee -a ${haplotype_id}.clean_stdout.log && CLEANED=true && break
    fi
  done

  # Fallback if no official script found
  if [ "\${CLEANED}" = "false" ]; then
    echo "No official FCS cleanup script found, using manual parser" | tee -a ${haplotype_id}.clean_stdout.log
    
    python3 << 'PYEOF'
import sys

# Simple FASTA parser (no BioPython dependency)
def read_fasta(filename):
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Get ID without description
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def write_fasta(sequences, filename):
    with open(filename, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f'>{seq_id}\\n')
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\\n')

# Parse action report to find sequences to exclude
exclude_seqs = set()
with open("${action_report}", 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        
        parts = line.strip().split('\\t')
        if len(parts) < 5:
            continue
        
        seq_id = parts[0]
        action = parts[3].upper()
        
        # EXCLUDE or TRIM actions mean contamination
        if 'EXCLUDE' in action or 'REMOVE' in action:
            exclude_seqs.add(seq_id)
            print(f"Excluding: {seq_id} (action: {action})", file=sys.stderr)

# Read and split sequences
all_seqs = read_fasta("${assembly_fa}")
clean_seqs = {k: v for k, v in all_seqs.items() if k not in exclude_seqs}
contam_seqs = {k: v for k, v in all_seqs.items() if k in exclude_seqs}

# Write outputs
write_fasta(clean_seqs, "${haplotype_id}.decontaminated.fasta")
write_fasta(contam_seqs, "${haplotype_id}.contaminants.fasta")

print(f"Clean sequences: {len(clean_seqs)}", file=sys.stderr)
print(f"Contaminated sequences: {len(contam_seqs)}", file=sys.stderr)
PYEOF

    echo "Manual cleanup completed" | tee -a ${haplotype_id}.clean_stdout.log
    echo "Clean sequences: \$(grep -c '^>' ${haplotype_id}.decontaminated.fasta || echo 0)" | tee -a ${haplotype_id}.clean_stdout.log
    echo "Contaminated sequences: \$(grep -c '^>' ${haplotype_id}.contaminants.fasta || echo 0)" | tee -a ${haplotype_id}.clean_stdout.log
  fi
  """
}