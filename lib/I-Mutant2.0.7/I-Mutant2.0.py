##  Copyright (C) 2004  Emidio Capriotti
## 
##  This program and all program in this package are free software;
##  you can redistribute it and/or modify it under the terms of the
##  GNU General Public License as published by the Free Software 
##  Foundation; either version 2 of the License, or (at your option)
##  any later version.
## 
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA
## 
##  contacts: Emidio Capriotti
##	      c/o Prof. Rita Casadio		
##            e-mail: emidio@biocomp.unibo.it
##            Dept. of Biology
##            University of Bologna
##            via Irnerio 42
##            40126 Bologna
##            Italy
## 
##  Last Update September 2004	
##
#import psyco

lh='VLIMFWYGAPSTCHRKQEND'

prese='\n'+75*'*'+'\n**'+71*' '+'**\n'+\
      '**'+29*' '+ 'I-Mutant v2.0'+29*' '+\
      '**\n**'+9*' '+'Predictor of Protein Stability Changes upon Mutations'+\
      9*' '+'**\n**'+71*' '+'**\n'+75*'*'+'\n\n'

leg= '      WT: Aminoacid in Wild-Type Protein\n'+\
     '      NEW: New Aminoacid after Mutation\n'+\
     '      DDG: DG(NewProtein)-DG(WildType) in Kcal/mole\n'+\
     '           DDG>0: Decrease Stability\n'+\
     '           DDG<0: Increase Stability\n'+\
     '      RI: Reliability Index\n'+\
     '      T: Temperature in Celsius unit\n'+\
     '      pH: -log[H+]\n'+\
     '      RSA: Relative Solvent Accessible Area\n'

cite= '\n'+75*'*'+'\n'+'*'+73*' '+'*'+'\n'+\
      '* http://gpcr.biocomp.unibo.it/cgi/predictors/I-Mutant2.0/I-Mutant2.0.cgi'+\
      ' '+'*\n'+'*'+73*' '+'*'+'\n'+75*'*'+'\n'

if (__name__  == '__main__'):
    import sys,string
    vinp=sys.argv
    cmd=vinp[0]
    ##print prese
    progpath=cmd[:cmd.find('I-Mutant2.0.py')]
    try:
	PDBT=progpath+'Tools/'
    	svm3dpath=progpath+'SvmData/3D/'
	svmseqpath=progpath+'SvmData/Seq/'
    	sys.path.append(PDBT)
	sys.path.append(svm3dpath)
	sys.path.append(svmseqpath)
	from Mutant import Svm3DOutput, SvmSeqOutput
	#psyco.bind(Svm3DOutput)
	#psyco.bind(SvmSeqOutput)
    except:
    	PDBT=progpath+'Tools\\'
	svm3dpath=progpath+'SvmData\\3D\\'
	svmseqpath=progpath+'SvmData\\Seq\\'
    	sys.path.append(PDBT)
	sys.path.append(svm3dpath)
	sys.path.append(svmseqpath)
	from Mutant import Svm3DOutput, SvmSeqOutput
	#psyco.bind(Svm3DOutput)
        #psyco.bind(SvmSeqOutput) 
    if len(vinp)==9 and vinp[1][:4]=='-pdb':
		pdbfile=vinp[2]
    		dsspfile=vinp[3]
    		chain=vinp[4]
    		position=vinp[5]
		try:
			inres=string.upper(vinp[6])
			if inres!='all' and inres!='ALL':
		        	if lh.find(inres[0])!=-1: newres=inres[0]
			else:
				newres='all'
		except:
		        newres='all'
		try:
			pH=float(vinp[7])
			if pH<=0.0 or pH>=100.0: pH=7.0
		except:
			pH=7.0
		try:
    			T=float(vinp[8])
			if T<=0.0 or T>=100.0: T=25.0
		except:
			T=25.0
		if vinp[1]!='-pdbv':
			opt='sddg'
		else:
			opt='ddg'
		vout=Svm3DOutput(pdbfile,dsspfile,chain,position,newres,opt,pH,T)
		
    elif len(vinp)==7  and vinp[1][:4]=='-pdb':
		pdbfile=vinp[2]
                dsspfile=vinp[3]
                chain=vinp[4]
                position=vinp[5]
		try:
			inres=string.upper(vinp[6])
			if inres!='all' and inres!='ALL':
				if lh.find(inres[0])!=-1: newres=inres[0]
			else:
				newres='all'
		except:
		        newres='all'
		if vinp[1]!='-pdbv':
                        opt='sddg'
                else:
                        opt='ddg'
		vout=Svm3DOutput(pdbfile,dsspfile,chain,position,newres,opt)

    elif len(vinp)==7 and  vinp[1][:4]=='-seq':
		seqfile=vinp[2]
		position=vinp[3]
                try:
			inres=string.upper(vinp[4])
			if inres!='all' and inres!='ALL':
		        	if lh.find(inres[0])!=-1: newres=inres[0]
			else:
				newres='all'
		except:
		        newres='all'
		try:
			pH=float(vinp[5])
			if pH<=0.0 or pH>=100.0: pH=7.0
		except:
			pH=7.0
		try:
    			T=float(vinp[6])
			if T<=0.0 or T>=100.0: T=25.0
		except:
			T=25.0
		if vinp[1]!='-seqv':
                        opt='sddg'
                else:
                        opt='ddg'
		vout=SvmSeqOutput(seqfile,position,newres,opt,pH,T)

    elif len(vinp)==5 and  vinp[1][:4]=='-seq':
		seqfile=vinp[2]
		position=vinp[3]
		try:
			inres=string.upper(vinp[4])
			if inres!='all' and inres!='ALL':
				if lh.find(inres[0])!=-1: newres=inres[0]
			else:
				newres='all'
		except:
		        newres='all'
		if vinp[1]!='-seqv':
                        opt='sddg'
                else:
                        opt='ddg'
		vout=SvmSeqOutput(seqfile,position,newres,opt)

    else:
		print '\nPrediction Using Structural Information:' 
		print ' - Predict Stability Change Direction'
		print '   python I-Mutant2.0.py -pdb  pdbfile dsspfile chain position newres pH=7.0 Temperature=25.0 '
		print ' - Predict Stability Change Value'
		print '   python I-Mutant2.0.py -pdbv pdbfile dsspfile chain position newres pH=7.0 Temperature=25.0'
		print 'Prediction Using Sequence Information:'
		print ' - Predict Stability Change Direction'
		print '   python I-Mutant2.0.py -seq  seqfile position newres pH=7.0 Temperature=25.0'
		print ' - Predict Stability Change Value'
		print '   python I-Mutant2.0.py -seqv seqfile position newres pH=7.0 Temperature=25.0\n'
		vexit=''
		while string.upper(vexit)!='Y':
                        kpred=''
			opt=''
			while kpred!='SEQ' and kpred!='PDB':
				kpred=string.upper(raw_input('Predictor Input (PDB=Structure, SEQ=Sequence)='))
			while opt=='':
				opt=string.upper(raw_input('Stability Output (S=Sign, V=Value)='))
				print opt
			if kpred=='PDB':
				pdbfile=''
				dsspfile=''
				while  pdbfile=='':
					pdbfile=raw_input('PDB file name: ')
				while dsspfile=='':
					dsspfile=raw_input('DSSP file name: ')
				chain=raw_input('PDB chain [_]: ')
				if chain=='' or len(chain)>1: chain='_'
			if kpred=='SEQ':
				seqfile=''
				seqfile=raw_input('Sequence file name: ')
			if opt!='V':
				opt='sddg'
			else:
				opt='ddg'
			position=''
			while position=='':
				try:
					position=int(raw_input('PDB residue position: '))
					position=str(position)
				except:
					position=''
			newres=''
			while newres=='':
				try:
					inres=string.upper(raw_input('New residue: '))
					if inres!='all' and inres!='ALL':
						if lh.find(inres[0])!=-1: newres=inres[0]
					else:
						newres='all'
				except:
					newres='all'
			try:
				pH=float(raw_input('pH [7.0]: '))
				if pH<=0.0 or pH>=14.0:
					pH=float(raw_input('pH (0.0-14.0): '))
			except:
				pH=7.0
			try:
				T=float(raw_input('Temperature [25.0 Celsius]: '))
				if T<=0.0 or T>=100.0:
					T=float(raw_input('Temperature (0-100 Celsius): '))
			except:
				T=25.0
			
			if kpred=='PDB': vout=Svm3DOutput(pdbfile,dsspfile,chain,position,newres,opt,pH,T)
			if kpred=='SEQ': vout=SvmSeqOutput(seqfile,position,newres,opt,pH,T)
			vexit=raw_input('Exit [N/Y]: ') 
    ##print '\n'+leg
    ##print cite
