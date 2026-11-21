# -*- coding: utf-8 -*-
"""Make the donor-level test primary and document the E-value conversion."""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from docx import Document

F = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\01_Manuscript.docx"
d = Document(F)
done = []


def sub(old, new, tag):
    for p in d.paragraphs:
        if old in p.text:
            p.runs[0].text = p.text.replace(old, new, 1)
            for r in p.runs[1:]:
                r.text = ''
            print('  %-14s OK' % tag)
            done.append(tag)
            return True
    print('  %-14s MISSING' % tag)
    return False


# ---------- 1. Results: donor level becomes primary ----------
OLD = ('VDR detection was aggregated by donor and cell-state group to avoid pseudoreplication '
       '(Fig. 4A, B), and modelled as detected versus undetected counts with a quasibinomial '
       'framework including donor identity and sequencing-depth metrics. In this adjusted model, '
       'RegC/HomC states remained enriched for VDR detection (OR = 1.89, 95% CI 1.32-2.69, '
       'P = 0.010; Table 4 and S3 Table). A paired donor-level Wilcoxon comparison showed the same '
       'median direction but did not reach significance (median RegC/HomC 13.2% versus Other 9.8%; '
       'P = 0.407), indicating donor-level heterogeneity.')
NEW = ('VDR detection was aggregated by donor and cell-state group so that the donor, not the cell, is '
       'the unit of inference (Fig. 4A, B). Across the ten donors of GSE104782 the RegC/HomC excess was '
       'positive in six, with a mean within-donor difference of 1.4 percentage points (95% CI -2.1 to '
       '4.9; paired Wilcoxon P = 0.375), so this atlas does not establish a donor-level difference. In '
       'GSE324993 the same contrast was positive in all six donors (mean +3.2 percentage points; '
       'Wilcoxon P = 0.031), and in GSE152805 the adjusted model gave OR 1.50 (1.09-2.06). A '
       'quasibinomial model that pools cells within donor-state units, with donor identity and '
       'sequencing-depth metrics as covariates, gave OR 1.89 (1.32-2.69, P = 0.010) for GSE104782; '
       'because that model treats cells as replicate observations, we report it as a supporting '
       'estimate rather than as donor-level evidence (Table 4 and S3 Table).')
sub(OLD, NEW, 'donor-primary')

# ---------- 2. cross-atlas paragraph: state where the donor signal is ----------
sub('The direction was therefore consistent across the three atlases: VDR detection was higher in '
    'regulatory states in GSE104782 (adjusted OR 1.89, 95% CI 1.32-2.69), GSE152805 (1.50, 1.09-2.06) '
    'and GSE324993 (1.13, 1.07-1.19).',
    'The direction was therefore consistent across the three atlases: VDR detection was higher in '
    'regulatory states in GSE104782 (adjusted OR 1.89, 95% CI 1.32-2.69), GSE152805 (1.50, 1.09-2.06) '
    'and GSE324993 (1.13, 1.07-1.19). At the donor level, however, the difference reached significance '
    'only in GSE324993 (six of six donors; Wilcoxon P = 0.031) and not in GSE104782 (six of ten; '
    'P = 0.375), so the replication is consistent in direction but not uniformly established at the '
    'donor level.',
    'cross-atlas')

# ---------- 3. Figure 4 caption: name the right test ----------
sub('(B) Within-donor difference in detection fraction; all six donors are positive, and the formal '
    'test of this contrast is the donor-adjusted quasibinomial model reported in the text.',
    '(B) Within-donor difference in detection fraction; all six donors are positive, and the formal '
    'test of this contrast is the paired donor-level Wilcoxon test reported in the text.',
    'fig4-caption')

# ---------- 4. Methods: E-value conversion ----------
sub('Survey-weighted models for the population-based cohort were fitted with the R survey package',
    'E-values were computed on the risk-ratio scale after converting the odds ratio with the outcome '
    'prevalence in the reference group (RR = OR / [1 - p0 + p0 x OR], p0 = 0.367), following '
    'VanderWeele and Ding[35]. Survey-weighted models for the population-based cohort were fitted with '
    'the R survey package',
    'evalue-method')

d.save(F)
print('\nedits:', len(done))
Q = [p.text for p in Document(F).paragraphs]
s = next(i for i, t in enumerate(Q) if t.strip().lower() == 'introduction')
r = next(i for i, t in enumerate(Q) if t.strip() == 'References')
print('abstract', len(Q[7].split()), '| main', len(' '.join(Q[s:r]).split()))
