# a script to calculate the best hits from a reciprocal protein blast 
# will produce some nice figures and output a tab delimited list of matchinh gene ids

#PACKAGES
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from pathlib import Path
import argparse
from Bio.Blast.Applications import NcbiblastpCommandline
from matplotlib.colors import LogNorm


def run_reciprocal_BLAST(query_aa, subject_aa, fwd_out, rev_out):
    # Create BLAST command-lines for forward and reverse BLAST searches
    fwd_blastp = NcbiblastpCommandline(query=query_aa, subject=subject_aa, out=fwd_out,
                                    outfmt='6 qseqid sseqid pident qcovs qlen slen length bitscore evalue',
                                    max_target_seqs=1)
    rev_blastp = NcbiblastpCommandline(query=subject_aa, subject=query_aa, out=rev_out,
                                    outfmt='6 qseqid sseqid pident qcovs qlen slen length bitscore evalue',
                                    max_target_seqs=1)
    
    # Inspect command-lines
    print("\033[46m {}\033[0;0m".format('BLAST searches to be run: '))
    print("FORWARD: %s" % fwd_blastp)
    print()
    print("REVERSE: %s" % rev_blastp)
    print()
    
    # run BLAST searches
    print("\033[32m {}\033[0;0m".format("Running forward BLAST..."))
    fwd_stdout, fwd_stderr = fwd_blastp()
    print("\033[32m {}\033[0;0m".format("Running reverse BLAST..."))
    rev_stdout, rev_stderr = rev_blastp()
    print("\033[32m {}\033[0;0m".format("BLAST complete!"))

    return

def normalise_blast_results(fwd_out, rev_out):
    # Load the BLAST results into Pandas dataframes
    fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
    rev_results = pd.read_csv(rev_out, sep="\t", header=None)

    headers = ["query", "subject", "identity", "coverage",
           "qlength", "slength", "alength",
           "bitscore", "E-value"]
    fwd_results.columns = headers
    rev_results.columns = headers

    # Create a new column in both dataframes: normalised bitscore
    fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
    rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

    # Create query and subject coverage columns in both dataframes
    fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
    rev_results['qcov'] = rev_results.alength/rev_results.qlength
    fwd_results['scov'] = fwd_results.alength/fwd_results.slength
    rev_results['scov'] = rev_results.alength/rev_results.slength

    # Clip maximum coverage values at 1.0
    fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
    rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
    fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
    rev_results['scov'] = rev_results['scov'].clip(upper=1)

    return fwd_results, rev_results

def identify_reciprocal_best_matches(fwd_results, rev_results, outName, outdir):
    # Merge forward and reverse results
    rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
                    left_on='subject', right_on='query',
                    how='outer')
    # Discard rows that are not RBH
    rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]
    # Group duplicate RBH rows, taking the maximum value in each column
    rbbh = rbbh.groupby(['query_x', 'subject_x']).max()
    path = os.path.join(outdir,outName+"_rbbh.tsv")
    rbbh.to_csv(path, index=None, header=True)
    return rbbh

def make_plots(fwd_results, rev_results, rbbh, outName,outdir):
    # Plot histograms
    f, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True)
    sns.despine(left=True)
    sns.distplot(fwd_results.norm_bitscore, color="b", ax=axes[0], axlabel="forward normalised bitscores")
    sns.distplot(rev_results.norm_bitscore, color="g", ax=axes[1], axlabel="reverse normalised bitscores")
    path = os.path.join(outdir,outName+'_norm_bitscores_histograms.png')
    plt.savefig(path, dpi=150)
    # Plot 2D density histograms
    '''
    Hits 'on the diagonal' have approximately the same coverage in query \
    and subject sequences: these are likely diverged proteins

    Hits 'off the diagonal' have more coverage in either the query or \
    subject sequence: these may be single-domain matches, poor alignments, \
    or some other result that is unlikely to be a very good match.
    '''
    # Calculate 2D density histograms for counts of matches at several coverage levels
    (Hfwd, xedgesf, yedgesf) = np.histogram2d(fwd_results.qcov, fwd_results.scov, bins=20)
    (Hrev, xedgesr, yedgesr) = np.histogram2d(rev_results.qcov, rev_results.scov, bins=20)
    # Create a 1x2 figure array
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharex=True, sharey=True)
    # Plot histogram for forward matches
    im = axes[0].imshow(Hfwd, cmap=plt.cm.Blues, norm=LogNorm(),
                        extent=[xedgesf[0], xedgesf[-1], yedgesf[0], yedgesf[-1]],
                        origin='lower', aspect=1)
    axes[0].set_title("Forward")
    axes[0].set_xlabel("query")
    axes[0].set_ylabel("subject")
    # Plot histogram for reverse matches
    im = axes[1].imshow(Hrev, cmap=plt.cm.Blues, norm=LogNorm(),
                        extent=[xedgesr[0], xedgesr[-1], yedgesr[0], yedgesr[-1]],
                        origin='lower', aspect=1)
    axes[1].set_title("Reverse")
    axes[1].set_xlabel("query")
    axes[1].set_ylabel("subject")
    # Add colourbars
    fig.colorbar(im, ax=axes[0])
    fig.colorbar(im, ax=axes[1])
    path = os.path.join(outdir,outName+"_norm_bitscores_2D_density_histogram.png")
    plt.savefig(path, dpi=150)

    # Plot distribution of RBH bitscores
    fig, axes = plt.subplots(1, 1, figsize=(6, 6), sharex=True, sharey=True)
    sns.distplot(rbbh.norm_bitscore, color="b", axlabel="RBH normalised bitscores")
    path = os.path.join(outdir,outName+"_RBH_norm_bitscores_histogram.png")
    plt.savefig(path, dpi=150)

    # Plot 2D density histograms for rbbh
    # Calculate 2D density histograms for counts of matches at several coverage levels
    (H, xedges, yedges) = np.histogram2d(rbbh.qcov, rbbh.scov, bins=20)

    # Create a 1x2 figure array
    fig, ax = plt.subplots(1, 1, figsize=(6, 6), sharex=True, sharey=True)

    # Plot histogram for RBBH
    im = ax.imshow(H, cmap=plt.cm.Blues, norm=LogNorm(),
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                    origin='lower', aspect=1)
    ax.set_title("RBBH")
    ax.set_xlabel("query")
    ax.set_ylabel("subject")

    # Add colourbar
    fig.colorbar(im, ax=ax)
    path = os.path.join(outdir,outName+"_RBH_norm_bitscores_2D_density_histogram.png")
    plt.savefig(path, dpi=150)
    return

def main():
    '''
    executes the arguments and functions
    adapted from: https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#imports
    '''
    parser = argparse.ArgumentParser(
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''\
Find matching gene models that might share gene IDs.
-------------------------------------------------------------

     ''',
     epilog="Adapted from: https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#imports by Helen Rebecca Davison") 
    parser.add_argument('-s1','--query_aa', \
                    help="one set of aa sequence",
                    required=True)
    parser.add_argument('-s2','--subject_aa', \
                    help="second set of aa sequence",
                    required=True)
    parser.add_argument('-o','--output', \
                    help="Name for outputs",
                    required=True)
    args = parser.parse_args()

    # inputs and names
    query_aa = args.query_aa
    subject_aa = args.subject_aa
    outName = args.output

    # define outputs
    outdir = os.path.join(os.getcwd(), outName+"-RBH-results")
    os.makedirs(outdir, exist_ok=True)
    fwd_out = os.path.join(outdir, outName+'-fwd-results.tab')
    rev_out = os.path.join(outdir, outName+'-rev-results.tab')

    # run the functions
    run_reciprocal_BLAST(query_aa, subject_aa, fwd_out, rev_out)
    fwd_results, rev_results = normalise_blast_results(fwd_out, rev_out)
    rbbh = identify_reciprocal_best_matches(fwd_results, rev_results, outName, outdir)
    make_plots(fwd_results, rev_results, rbbh, outName, outdir)
    return
main()
