# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"

import unittest, os, tempfile, shutil, subprocess, argparse
import taxon_filter, util.file, tools.last, tools.bmtagger, tools.blast
from test import assert_equal_contents, TestCaseWithTmp


class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in taxon_filter.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestTrimmomatic(TestCaseWithTmp) :

    def test_trimmomatic(self) :
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        pairedOutFastq1 = util.file.mkstempfname()
        pairedOutFastq2 = util.file.mkstempfname()
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        parser = taxon_filter.parser_trim_trimmomatic(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, pairedOutFastq1,
            pairedOutFastq2, clipFasta])
        args.func_main(args)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        assert_equal_contents(self, pairedOutFastq1, expected1Fastq)
        assert_equal_contents(self, pairedOutFastq2, expected2Fastq)

class TestFilterLastal(TestCaseWithTmp) :

    def test_filter_lastal(self) :
        # Create refDbs
        commonInputDir = util.file.get_test_input_path()
        myInputDir = util.file.get_test_input_path(self)
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')
        dbsDir = tempfile.mkdtemp()
        refDbs = os.path.join(dbsDir, 'ebola')
        lastdbPath = tools.last.Lastdb().install_and_get_path()
        subprocess.check_call([lastdbPath, refDbs, refFasta])

        # Call main_filter_lastal
        inFastq = os.path.join( myInputDir, 'in.fastq')
        outFastq = util.file.mkstempfname('.fastq')
        args = taxon_filter.parser_filter_lastal(argparse.ArgumentParser()).parse_args([
            inFastq, refDbs, outFastq])
        args.func_main(args)

        # Check that results match expected
        expectedFastq = os.path.join(myInputDir, 'expected.fastq')
        assert_equal_contents(self, outFastq, expectedFastq)

class TestFilterLastalBam(TestCaseWithTmp) :

    def test_filter_lastal(self) :
        # Create refDbs
        commonInputDir = util.file.get_test_input_path()
        myInputDir = util.file.get_test_input_path(self)
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')
        dbsDir = tempfile.mkdtemp()
        refDbs = os.path.join(dbsDir, 'ebola')
        lastdbPath = tools.last.Lastdb().install_and_get_path()
        subprocess.check_call([lastdbPath, refDbs, refFasta])

        # create bam from test fastq's
        testInputBam = util.file.mkstempfname('.bam')
        tools.picard.FastqToSamTool().execute(os.path.join(myInputDir, 'in.fastq'),
                                              os.path.join(myInputDir, 'in2.fastq'),
                                              "TestFilterLastalBam", testInputBam,
                                              ["V=Standard"])

        # Call parser_filter_lastal_bam
        testOutputBam = util.file.mkstempfname('.bam')
        args = taxon_filter.parser_filter_lastal_bam(argparse.ArgumentParser()).parse_args([
            testInputBam, refDbs, testOutputBam])
        args.func_main(args)

        # convert bam back to fastq's
        testOutputFastq1 = util.file.mkstempfname('.fastq')
        testOutputFastq2 = util.file.mkstempfname('.fastq')
        tools.picard.SamToFastqTool().execute(testOutputBam, testOutputFastq1, testOutputFastq2)

        # Check that results match expected
        assert_equal_contents(self, testOutputFastq1, os.path.join(myInputDir, 'expected.fastq'))
        assert_equal_contents(self, testOutputFastq2, os.path.join(myInputDir, 'expected2.fastq'))

class TestBmtagger(TestCaseWithTmp) :
    """
    How test data was created:
      humanChr1Subset.fa has 200 bases from human chr1
      humanChr9Subset.fa has 200 bases from human chr9
      bmtool -d humanChr1Subset.fa -o humanChr1Subset.bitmask -w 8
      bmtool -d humanChr9Subset.fa -o humanChr9Subset.bitmask -w 8
      in[12].fastq "reads" are from humanChr[19]Subset.fa and ebola genome,
          with arbitrary quality scores.
    """
    def setUp(self) :
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        srprismPath = tools.bmtagger.SrprismTool().install_and_get_path()
        for db in ['humanChr1Subset', 'humanChr9Subset'] :
            # .map file is > 100M, so recreate instead of copying
            dbfa = os.path.join(myInputDir, db + '.fa')
            dbsrprism = os.path.join(self.tempDir, db + '.srprism')
            subprocess.check_call([srprismPath, 'mkindex', '-i', dbfa, '-o', dbsrprism])
            # .bitmask and .srprism.* files must be in same dir, so copy
            shutil.copy(os.path.join(myInputDir, db + '.bitmask'), self.tempDir)

    def test_partition_bmtagger(self) :
        outMatch   = [os.path.join(self.tempDir,   'outMatch.{}.fastq'.format(n))
                      for n in '12']
        outNoMatch = [os.path.join(self.tempDir, 'outNoMatch.{}.fastq'.format(n))
                      for n in '12']
        myInputDir = util.file.get_test_input_path(self)
        args = taxon_filter.parser_partition_bmtagger(argparse.ArgumentParser()).parse_args(
            [os.path.join(myInputDir, 'in1.fastq'),
             os.path.join(myInputDir, 'in2.fastq'),
             os.path.join(self.tempDir, 'humanChr1Subset'),
             os.path.join(self.tempDir, 'humanChr9Subset'),
             '--outMatch', outMatch[0], outMatch[1],
             '--outNoMatch', outNoMatch[0], outNoMatch[1]])
        args.func_main(args)
            
        # Compare to expected
        for case in ['Match.1', 'Match.2', 'NoMatch.1', 'NoMatch.2'] :
            assert_equal_contents(self,
                os.path.join(self.tempDir, 'out' + case + '.fastq'),
                os.path.join(myInputDir, 'expected.' + case + '.fastq'))

    def test_deplete_bmtagger(self) :
        myInputDir = util.file.get_test_input_path(self)
        args = taxon_filter.parser_partition_bmtagger(argparse.ArgumentParser()).parse_args(
            [os.path.join(myInputDir, 'in1.fastq'),
             os.path.join(myInputDir, 'in2.fastq'),
             os.path.join(self.tempDir, 'humanChr1Subset'),
             os.path.join(self.tempDir, 'humanChr9Subset'),
             '--outNoMatch',
             os.path.join(self.tempDir, 'deplete.1.fastq'),
             os.path.join(self.tempDir, 'deplete.2.fastq')])
        args.func_main(args)
        
        # Compare to expected
        for case in ['1', '2'] :
            assert_equal_contents(self,
                os.path.join(self.tempDir, 'deplete.' + case + '.fastq'),
                os.path.join(myInputDir, 'expected.NoMatch.' + case + '.fastq'))

class TestDepleteBlastn(TestCaseWithTmp) :
    '''
    How test data was created:
      humanChr1Subset.fa has 200 bases from human chr1
      humanChr9Subset.fa has 200 bases from human chr9
      in.fastq "reads" are from humanChr[19]Subset.fa and ebola genome,
          with arbitrary quality scores.
    '''
    def test_deplete_blastn(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)

        # Make blast databases
        makeblastdbPath = tools.blast.MakeblastdbTool().install_and_get_path()
        dbnames = ['humanChr1Subset.fa', 'humanChr9Subset.fa']
        refDbs = []
        for dbname in dbnames :
            refDb = os.path.join(tempDir, dbname)
            os.symlink(os.path.join(myInputDir, dbname), refDb)
            refDbs.append(refDb)
            subprocess.check_call([
                makeblastdbPath, '-dbtype', 'nucl', '-in', refDb])

        # Run deplete_blastn
        outFile = os.path.join(tempDir, 'out.fastq')
        args = taxon_filter.parser_deplete_blastn(argparse.ArgumentParser()).parse_args(
            [os.path.join(myInputDir, 'in.fastq'),
             outFile,
             refDbs[0],
             refDbs[1]])
        args.func_main(args)

        # Compare to expected
        assert_equal_contents(self, outFile,
                              os.path.join(myInputDir, 'expected.fastq'))


if __name__ == '__main__':
    unittest.main()
