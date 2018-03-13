"""

pymzid - python mzIdentML Parser v.0.3.0.
pymzid reads in mzid files and creates flat summary tables.
Written by Edward Lau (lau1@stanford.edu) 2016-2018

Example:
    parse_mzid.py percolator.target.mzid --out=mzid.txt

"""

from time import time
from xml.dom import minidom
import pandas as pd
import os.path
import sys
import re
import numpy as np

class Mzid(object):
    """
    This class holds the loaded mzIdentml object and extracts peptides and protein labels.
    Goal is to get dataframe where identified peptide sequence/z combos are linked to their m/z, protein and spectrum ID
    To do so, first pull 5 tables from the mzID file: 1.PSM, 2.Peptide, 3.ProteinGroup, 4.PeptideEvidence, 5.DBSequence

    1. PSM contains spectrum_ID, mz, z passThreshold; references to Peptide and PeptideEvidence through Pep_ref & PE_ref
    2. Peptide contains the actual amino acid sequence
    3. ProteinGroup contains passThreshold, references to DBSequence and PeptideEvidence through DBS_ref and PE_Ref
    4. PeptideEvidence contains isDecoy, references to Peptide and DBSequence through Pep_ref and DBS_ref
    5. DBSequence contains the Uniprot accession number

    From these five tables, output a peptide-centric summary, and a protein-centric summary
    Peptide-centric Summary combines PSM and Peptide, contains Sequence, mz, z, spectrumID
    Protein-centric Summary reads all 5 tables, should contain Uniprot in addition to peptide-centric summary

    """

    def __init__(self, path):
        """
        :param path: path of the mzid file to be loaded, e.g., "~/Desktop/example.mzid"


        """
        self.path = os.path.join(path)
        self.root = self.parse_file()

        self.psm_df = self.read_psm()
        self.peptide_df = self.read_peptide()
        self.protein_df = self.read_protein()
        self.pe_df = self.read_pe()
        self.dbs_df = self.read_dbs()

        self.pep_summary_df = pd.DataFrame()
        self.pro_summary_df = pd.DataFrame()

        self.filtered_protein_df = pd.DataFrame()
        self.filtered_pep_summary_df = pd.DataFrame()

    def parse_file(self):
        """
        Get the mzid file xml root
        NB: Move time tracker to main loop when made into object oriented.

        :return: True
        """

        print('Reading mzID file as document object model...')
        t1 = time()
        xmldoc = minidom.parse(self.path).childNodes[0]
        # xmldoc = minidom.parse("small_test/percolator.target.mzid").childNodes[0]
        # xmldoc = minidom.parse("external_mzid_test/mzidentml-example.mzid").childNodes[0]
        # xmldoc = minidom.parse("external_mzid_test/BSA1_msgfplus_v2016_09_16.mzid").childNodes[0]

        t2 = time()
        print('Done. Processing time: ' + str(round(t2 - t1, 2)) + ' seconds.')

        # Check if this is mzML 1.1 or above
        if xmldoc.attributes['version'].value < '1.1.0':
            sys.exit('mzML version is not 1.1.0 or above.')

        else:
            return xmldoc

    def read_psm(self):
        """
        Table #1 : All PSM Identified

        :return: psm_df     - A dataframe containing all peptide spectrum matches
        """

        allPsmList = []

        # Traverse to DataCollection > SpectrumIdentificationList
        # Note we jump past AnalysisData, assuming there is only one such section
        DataCollection = self.root.getElementsByTagName('DataCollection')[0]
        SpectrumIdentificationList = DataCollection.getElementsByTagName('SpectrumIdentificationList')[0]

        # How many SpectrumIdentificationResult entries are in each SpectrumIdentificationList ?
        # Traverse down to SpectrumIdentificationResult. Each SIR is linked to one spectrum ID but possibly multiple SII
        # Each SII is a PSM.
        SpectrumIdentificationResult = SpectrumIdentificationList.getElementsByTagName('SpectrumIdentificationResult')

        for i in range(0, len(SpectrumIdentificationResult)):

            # Print progress every 1000 records
            if i % 1000 == 0:
                print('Processing ' + str(i) + ' of ' + str(len(SpectrumIdentificationResult)) + ' records (PSM).')

            sir_id = SpectrumIdentificationResult[i].attributes['id'].value
            spectrum_id = SpectrumIdentificationResult[i].attributes['spectrumID'].value.split('-')[0]

            # Traverse down to the SII entries inside each SIR
            SpectrumIdentificationItem = SpectrumIdentificationResult[i].getElementsByTagName(
                'SpectrumIdentificationItem')

            for j in range(0, len(SpectrumIdentificationItem)):

                psmList = []

                sii_id = SpectrumIdentificationItem[j].attributes['id'].value
                z = SpectrumIdentificationItem[j].attributes['chargeState'].value
                mz = SpectrumIdentificationItem[j].attributes['experimentalMassToCharge'].value
                calc_mz = SpectrumIdentificationItem[j].attributes['calculatedMassToCharge'].value
                pep_id = SpectrumIdentificationItem[j].attributes['peptide_ref'].value
                pass_threshold = SpectrumIdentificationItem[j].attributes['passThreshold'].value

                # Get the Peptide Evidence Ref (need this to link to other tables later)
                try:
                    PeptideEvidenceRef = SpectrumIdentificationItem[j].getElementsByTagName('PeptideEvidenceRef')[0]
                    pe_id = PeptideEvidenceRef.attributes['peptideEvidence_ref'].value
                except IndexError:
                    pe_id = "no_id"

                # Put SIR ID, spectrum_ID, SII ID, z, m/z, theoretical m/z, Peptide ID, passThreshold, PE ID into list
                psmList.extend([sir_id, spectrum_id, sii_id, z, mz, calc_mz, pep_id, pass_threshold, pe_id])

                # Get the list of all cvParams from each SII (e.g.,percolator scores)
                cvParams = SpectrumIdentificationItem[j].getElementsByTagName('cvParam')

                for cvParam in cvParams:

                    # Restrict to Top level cvParams only (no fragmentation stuff)
                    if cvParam.parentNode != SpectrumIdentificationItem[j]:
                        continue

                    acc = cvParam.attributes['accession'].value
                    name = cvParam.attributes['name'].value
                    try:
                        value = cvParam.attributes['value'].value
                    except KeyError:
                        value = 0
                    cvParamList = [acc, name, value]
                    psmList.append(cvParamList)

                    # peptideList.append(acc)
                    # peptideList.append(name)
                    # peptideList.append(value)

                    # Add each PSM's psmList of info into the list of lists
                    allPsmList.append(psmList)



        # Convert into pandas dataframe
        psm_df = pd.DataFrame(allPsmList)
        psm_df.to_csv("all_psms.txt", sep='\t')
        print(psm_df)
        #
        # Rename the columns in the PSM table
        #

        # Get column names
        colnames = psm_df.columns.tolist()

        # Create a dictionary that map old names to new names
        # The first three columns are ID, Sequence
        newcol = {0: 'sir_id',
                  1: 'spectrum_id',
                  2: 'sii_id',
                  3: 'z',
                  4: 'mz',
                  5: 'calc_mz',
                  6: 'pep_id',
                  7: 'pep_pass_threshold',
                  8: 'pe_id',
                  }

        # Getting the name of the cvParam from the first filled row of the ith column of the dataframe
        for i in range(9, len(colnames)):
            for j in range(0, len(psm_df)):
                try:
                    newcol[i] = re.sub(' |:', '_', psm_df[i][j][1])
                    break
                except TypeError:
                    continue

            temp = []

            # Loop over each line and add values
            # This is rewritten to use [row][column] index instead of just adding the entire column
            # e.g., like for accession, name, value in psm[i] temp.append(value)
            # In order to account for some search engines having certain cvParam flags in only some rows.
            for j in range(0, len(psm_df)):
                if psm_df[i][j] is None:
                    temp.append(-1)
                else:
                    temp.append(psm_df[i][j][2])

            psm_df[i] = temp

        psm_df.rename(columns=newcol, inplace=True)

        return psm_df

    def read_peptide(self):
        """
        Table # 2: Peptide Sequences Identified

        Get the values inside Peptide. This is needed because the SII (PSM) refers only to Peptide_ref
        Peptide_ref needs to be looked up to get the sequence

        :return:
        """

        ## Traverse down to Sequence Collection > Peptide
        SequenceCollection = self.root.getElementsByTagName('SequenceCollection')[0]
        Peptide = SequenceCollection.getElementsByTagName('Peptide')

        allPeptideList = []

        for i in range(0, len(Peptide)):

            peptideList = []

            # Print progress every 1000 records
            if i % 1000 == 0:
                print('Processing ' + str(i) + ' of ' + str(len(Peptide)) + ' records (Peptides).')

            pep_id = Peptide[i].attributes['id'].value

            PeptideSequence = Peptide[i].getElementsByTagName('PeptideSequence')[0]
            seq = PeptideSequence.childNodes[0].nodeValue

            # Put Peptide ID and peptide sequence into peptideList
            peptideList.extend([pep_id, seq])

            # Get the cvParams from each Peptide (percolator scores)
            cvParams = Peptide[i].getElementsByTagName('cvParam')

            for cvParam in cvParams:

                # Restrict to Top level cvParams only (no fragmentation stuff)
                if cvParam.parentNode != Peptide[i]:
                    continue

                acc = cvParam.attributes['accession'].value
                name = cvParam.attributes['name'].value

                try:
                    value = cvParam.attributes['value'].value
                except KeyError:
                    value = 0

                cvParamList = [acc, name, value]
                peptideList.append(cvParamList)

            allPeptideList.append(peptideList)

        peptide_df = pd.DataFrame(allPeptideList)

        #
        # Rename the columns in the peptide table
        #

        # Get column names
        colnames = peptide_df.columns.tolist()

        # Create a dictionary that map old names to new names
        # The first three columns are ID, Sequence
        newcol = {0: 'pep_id',
                  1: 'seq',
                  }

        # Getting the name of the cvParam from the first row of the ith column of the dataframe
        for i in range(2, len(colnames)):
            newcol[i] = re.sub(' |:', '_', peptide_df[i][0][1])
            temp =[]
            for acc, name, value in peptide_df[i]:
                temp.append(value)
            peptide_df[i] = temp

        peptide_df.rename(columns=newcol, inplace=True)

        return peptide_df

    def read_protein(self):
        """
         Retrieve ProteinAmbiguityGroup for all identified proteins from the mzIdentML file.

        :return:
        """

        # Traverse to DataCollection > ProteinDetectionList
        DataCollection = self.root.getElementsByTagName('DataCollection')[0]
        ProteinDetectionList = DataCollection.getElementsByTagName('ProteinDetectionList')


        # If there is no ProteinDetectionList in the mzid (apparently seen in some MSGF results)
        # Return an empty protein_df
        if not ProteinDetectionList:
            print('No Protein Detection List found.')
            protein_df = pd.DataFrame(columns=['pag_id',
                                               'pdh_ud',
                                               'dbs_id',
                                               'prot_pass_threshold',
                                               'pe_id',
                                               'sii_id'])
            return protein_df

        # Traverse down to each ProteinAmbiguityGroup from each ProteinDetectionList - There are multiple PAG per PDL
        ProteinAmbiguityGroup = ProteinDetectionList[0].getElementsByTagName('ProteinAmbiguityGroup')

        allProteinList = []

        for i in range(0, len(ProteinAmbiguityGroup)):

            # Print progress every 1000 records
            if i % 1000 == 0:
                print('Processing ' + str(i) + ' of ' + str(len(ProteinAmbiguityGroup)) + ' records (Proteins).')

            pag_id = ProteinAmbiguityGroup[i].attributes['id'].value

            # Traverse down to each Protein Detection Hypothesis from each Protein Ambiguity Group.
            # There may be >1 PDH per PAG in some files, but in the test files I have done.
            # (From Comet and Tide, but only single fraction) there is only 1 PDH per PAG.
            # It is possible that for multi-fraction runs there can be multiple PDH per PAG.
            ProteinDetectionHypothesis = ProteinAmbiguityGroup[i].getElementsByTagName('ProteinDetectionHypothesis')[0]

            # Grab the ProteinDetectionHypothesis ID, DBSequence Reference, and Threshold flag
            pdh_id = ProteinDetectionHypothesis.attributes['id'].value
            dbs_id = ProteinDetectionHypothesis.attributes['dBSequence_ref'].value
            pass_threshold = ProteinDetectionHypothesis.attributes['passThreshold'].value


            # Get the cvParams from each ProteinDetectionHypothesis (e.g., percolator scores)
            cvParams = ProteinDetectionHypothesis.getElementsByTagName('cvParam')

            # Traverse down to each Peptide Hypothesis from each Protein Detection Hypothesis

            PeptideHypothesis = ProteinDetectionHypothesis.getElementsByTagName('PeptideHypothesis')

            for j in range(0, len(PeptideHypothesis)):
                pe_id = PeptideHypothesis[j].attributes['peptideEvidence_ref'].value

                # Traverse down to each SII Ref from each PeptideHypothesis
                SpectrumIdentificationItemRef = PeptideHypothesis[j].getElementsByTagName('SpectrumIdentificationItemRef')

                for k in range(0, len(SpectrumIdentificationItemRef)):
                    sii_id = SpectrumIdentificationItemRef[k].attributes['spectrumIdentificationItem_ref'].value

                    # Put pag_id, pdh_id, dbs_ref, pass_threshold, pe_ref, and sii_ref into the proteinList of a particular PAG
                    proteinList = []
                    proteinList.extend([pag_id, pdh_id, dbs_id, pass_threshold, pe_id, sii_id])

                    for cvParam in cvParams:

                        # Restrict to Top level cvParams only (no fragmentation stuff)
                        if cvParam.parentNode != ProteinDetectionHypothesis:
                            continue

                        acc = cvParam.attributes['accession'].value
                        name = cvParam.attributes['name'].value

                        try:
                            value = cvParam.attributes['value'].value
                        except KeyError:
                            value = 0

                        cvParamList = [acc, name, value]
                        proteinList.append(cvParamList)

                    allProteinList.append(proteinList)

        protein_df = pd.DataFrame(allProteinList)

        #
        # Rename the columns in the peptide table
        #

        # Get column names
        colnames = protein_df.columns.tolist()

        # Create a dictionary that map old names to new names
        # The first three columns are ID, Sequence
        newcol = {0: 'pag_id',
                  1: 'pdh_id',
                  2: 'dbs_id',
                  3: 'prot_pass_threshold',
                  4: 'pe_id',
                  5: 'sii_id',
                  }

        # Getting the name of the cvParam from the first row of the ith column of the dataframe
        for i in range(6, len(colnames)):
            newcol[i] = re.sub(' |:', '_', protein_df[i][0][1])
            temp = []
            for acc, name, value in protein_df[i]:
                temp.append(value)
            protein_df[i] = temp


        protein_df.rename(columns=newcol, inplace=True)

        return protein_df


    def read_pe(self):
        # Table #4. PeptideEvidence
        # This is read to see whether sequence isDecoy, and to link each PAG's PeptideHypothesis to a Peptide and a DBSequence


        ## Traverse down to Sequence Collection > DBSequence
        SequenceCollection = self.root.getElementsByTagName('SequenceCollection')[0]
        PeptideEvidence = SequenceCollection.getElementsByTagName('PeptideEvidence')

        allPEList = []

        for i in range(0, len(PeptideEvidence)):

            peList = []

            # Print progress every 1000 records
            if i % 1000 == 0:
                print('Processing ' + str(i) + ' of ' + str(len(PeptideEvidence)) + ' records (PeptideEvidence).')

            pe_id = PeptideEvidence[i].attributes['id'].value
            pep_id = PeptideEvidence[i].attributes['peptide_ref'].value
            dbs_id = PeptideEvidence[i].attributes['dBSequence_ref'].value
            start = PeptideEvidence[i].attributes['start'].value
            end = PeptideEvidence[i].attributes['end'].value
            is_decoy = PeptideEvidence[i].attributes['isDecoy'].value

            # Add PeptideEvidence ID and reference to Peptide, reference DBSequence, and isDecoy flag to peList
            peList.extend([pe_id, pep_id, dbs_id, start, end, is_decoy])

            allPEList.append(peList)

        pe_df = pd.DataFrame(allPEList)

        #
        # Rename the columns in the peptide table
        #

        # Get column names
        colnames = pe_df.columns.tolist()

        # Create a dictionary that map old names to new names
        # The first three columns are ID, Sequence
        newcol = {0: 'pe_id',
                  1: 'pep_id',
                  2: 'dbs_id',
                  3: 'start',
                  4: 'end',
                  5: 'is_decoy',
                  }

        pe_df.rename(columns=newcol, inplace=True)

        return pe_df


    def read_dbs(self):
        #
        # Table #5: DBSequences
        # This is read to get the Uniprot accession of each identified peptide
        #


        ## Traverse down to Sequence Collection > DBSequence
        SequenceCollection = self.root.getElementsByTagName('SequenceCollection')[0]
        DBSequence = SequenceCollection.getElementsByTagName('DBSequence')

        allDBSList = []

        for i in range(0, len(DBSequence)):

            dbsList = []

            # Print progress every 1000 records
            if i % 1000 == 0:
                print('Processing ' + str(i) + ' of ' + str(len(DBSequence)) + ' records (DatabaseSequences).')

            dbs_id = DBSequence[i].attributes['id'].value

            try:
                length = DBSequence[i].attributes['length'].value
            except KeyError:
                length = -1

            acc = DBSequence[i].attributes['accession'].value

            # Add DBSequence ID, sequence length, and Uniprot accession to the dbsList
            dbsList.extend([dbs_id, length, acc])

            allDBSList.append(dbsList)

        dbs_df = pd.DataFrame(allDBSList)

        #
        # Rename the columns in the peptide table
        #

        # Get column names
        colnames = dbs_df.columns.tolist()

        # Create a dictionary that map old names to new names
        # The first three columns are ID, Sequence
        newcol = {0: 'dbs_id',
                  1: 'length',
                  2: 'acc',
                  }

        dbs_df.rename(columns=newcol, inplace=True)

        return dbs_df

    def make_peptide_summary(self, take_psm_df_cvParams=True):
        """
        Take the five Mzid data frames, and return a flat table after filtering for Q values

        It first combines PSM and Peptide into a peptide-centric summary, which contains:
        Peptide ID, sequence, ... (cvParams in peptide_df) ... , SIR_id, spectrum ID, SII_Id
        z, m/z, calculated m/z, pep_pass_threshold, pe_id, ... (optionally, cvParams in psm_df)

        2016-03-06 EL
        Unfortunately the ugly "take_psm_df_cvParams" flag is needed here
        because some workflows like Percolator has cvParams with identical names for different purposes
        in both the PSM section and the Peptide section.

        Percolator outputs "percolator:score",
        "percolator:Q value", and "percolator:PEP" for both Spectral Identification Item and for Peptide
        In Riana we only take the peptide level Q value but not PSM-level q value, hence if the
        take_psm_df_cvParams flag is set to False we will take only the standard columns [0:8]
        from the psm_df (without the cvParams-derived columns)
        in order to avoid conflicts in merging.

        We may have to create other flags (such as to only take certain standard
        columns from peptide_df and other dataframes) to account for other cvParams conflict.

        :param take_psm_df_cvParams: Whether to keep the cvParams columns from the PSM dataframe
        :return: True
        """

        # Here we assuming the identification scores are already sorted from top to bottom
        # #and take only one psm_id for each pep_id.
        self.psm_df = self.psm_df.groupby('pep_id', sort=False).first().reset_index()

        if take_psm_df_cvParams:
            self.pep_summary_df = pd.merge(self.peptide_df, self.psm_df, how='left')
        else:
            self.pep_summary_df = pd.merge(self.peptide_df, self.psm_df[self.psm_df.columns[0:8]], how='left')

        return True

    def filter_peptide_summary(self, lysine_filter=0, protein_q=1e-2, peptide_q=1e-2, unique_only=False, require_protein_id=False):
        """
        The peptide-centric summary is then fitered by:
        - peptides that belong to any protein identified at a protein Q value
        - peptides below a certain peptide Q value
        - peptides containing certain number of lysines

        # 2017-04-06 Note the require_protein_id flag doesn't work for MSGF+ test files at the moment
        # because it seems the MSGF+ mzID files have no ProteinDetectioNList fields but instead
        # store the protein accessions inside <DBSequence>. Turn the flag off when doing MSGF+.

        :param lysine_filter: Lysine filter from command line argument
        :param protein_q: Protein-level Q value from command line argument
        :param peptide_q: Peptide-level Q value from command line argument
        :param unique_only: Only doing unique peptides
        :param require_protein_id: Require protein IDs (to filter out some mascot rows with no protein fields)
        :return: True
        """

        #
        # Filter peptides by Protein Q value.
        #
        try:
            self.filtered_protein_df = self.protein_df.loc[lambda x: x.percolator_Q_value.astype(float) < protein_q, :]
            #self.filtered_protein_df = self.filtered_protein_df.reset_index()
        except:
            print('No filtering by protein Q value done.')
            self.filtered_protein_df = self.protein_df


        #
        # Filter peptides by peptide-level Q-value filter
        #
        try:
            self.pep_summary_df = self.pep_summary_df.loc[lambda x: x.percolator_Q_value.astype(float) < peptide_q, :]
            self.pep_summary_df = self.pep_summary_df.reset_index(drop=True)

        except:
            pass

        #
        # Filter peptides by optional lysine filter
        #

        # If the lysine filter is on and the peptide has no or more than one lysine, skup
        if lysine_filter == 1:
            try:
                self.pep_summary_df = self.pep_summary_df.loc[lambda x: x.seq.apply(lambda y: y.count('K')) == 1, :]
                self.pep_summary_df = self.pep_summary_df.reset_index(drop=True)

            except:
                pass

        if lysine_filter == 2:
            try:
                self.pep_summary_df = self.pep_summary_df.loc[lambda x: x.seq.apply(lambda y: y.count('K')) > 0, :]
                self.pep_summary_df = self.pep_summary_df.reset_index(drop=True)

            except:
               pass

        # 2013-04-06 Again I forgot why we only chose the first five columns here. Resetting to all columns for now.

        if require_protein_id:
            self.filtered_pep_summary_df = pd.merge(self.pep_summary_df,
                                                    self.filtered_protein_df[self.filtered_protein_df.columns[[0, 5]]],
                                                    how='inner') # NB: 2017-04-05 might need 'inner' here for riana to work
        elif not require_protein_id:
            self.filtered_pep_summary_df = pd.merge(self.pep_summary_df,
                                                    self.filtered_protein_df[self.filtered_protein_df.columns[[0, 5]]],
                                                    how='left')  # NB: 2017-04-05 might need 'inner' here for riana to work



        #
        # Get the protein Uniprot accession via PE and then via DBS
        #
        self.filtered_pep_summary_df = pd.merge(self.filtered_pep_summary_df, self.pe_df, how='left')
        self.filtered_pep_summary_df = pd.merge(self.filtered_pep_summary_df, self.dbs_df, how='left')
        self.filtered_pep_summary_df = self.filtered_pep_summary_df.reset_index(drop=True)

        # Get only the peptides associated with one and only one proteins
        if unique_only:

            try:
                self.filtered_pep_summary_df = self.filtered_pep_summary_df.groupby('seq').filter(lambda x: (len(set(x['acc'])) == 1))
                self.filtered_pep_summary_df = self.filtered_pep_summary_df.reset_index(drop=True)

            except:
                pass

        self.filtered_pep_summary_df = self.filtered_pep_summary_df.reset_index(drop=True)

        return True

def read(args):
    """
    Parse
    :param args: Arguments from command line
    :return:

    """

    # Handle command line arguments
    mzid_loc = args.mzid
    out_loc = args.out
    prot_only = args.id

    try:
        mzid = Mzid(mzid_loc)
        mzid.make_peptide_summary(take_psm_df_cvParams=True)
        mzid.filter_peptide_summary(lysine_filter=0, protein_q=1, peptide_q=1, unique_only=False, require_protein_id=prot_only)
        mzid.filtered_pep_summary_df.to_csv(out_loc, sep='\t')


    except OSError as e:
        sys.exit('Failed to load mzid file. ' + str(e.errno))

    return sys.exit(os.EX_OK)


#
# Code for running main with parsed arguments from command line
#


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='PyMzid v.0.3.0 reads protein identification mzID files')

    parser.add_argument('mzid', help='path to mzid file')
    parser.add_argument('-i' '--id', action='store_true',
                        help='only outputs rows associated with protein accession')
    parser.add_argument('-o', '--out', help='name of the output files [default: mzid.txt]',
                        default='mzid.txt')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose error messages')


    parser.set_defaults(func=read)

    # Print help message if no arguments are given
    import sys
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)
