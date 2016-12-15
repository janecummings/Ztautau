import ROOT
import streams
from Drawer import *
import os,sys
from array import array

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetFrameBorderSize(1)

int_lumi = 20243.6


def compare_SR_truth( mu, entry, hand, norm = False):


	rf = ROOT.TFile('~/panalysis/output/prod4/theory/%s_regions.root' % mu.name)

	cname = 'Reweight_Pythia8_eta_%s_%s_%s' % (mu.name, entry.name, hand)
	pname = 'Reweight_Pythia_eta_SR_%s' % (mu.name)
	if norm:
		cname += '_norm'
		pname += '_norm'


	c = ROOT.TCanvas(cname,cname,1500,1000)
	c.UseCurrentStyle()
	c.Divide(1,2)
	c.cd(1).SetPad(0.0,0.2,1.0,1.0)
	c.cd(2).SetPad(0.0,0.0,1.0,0.2)
	c.cd(1)

	#l = ROOT.TLegend(.68,.7,.85,.9)
	l = ROOT.TLegend(.65,.76,.9,.9)
	l.SetTextSize(0.035)
	l.SetBorderSize(0)
	l.SetFillColor(0)

	region = 'SR_OS'

	hists = []
	maxi = 0
	col = 0
	for i,group in enumerate(mu.MC[0:2]+mu.Alpgen):
		if not hand in group.name:
			continue

		h = rf.Get('h_%s_%s_%s' % (entry.name, group.name, region))
		h.Rebin(entry.rebin)
		if norm:
			h.Scale( 1./ h.Integral())

		maxi = max( maxi, h.GetMaximum())
		h.SetLineColor( col + 1)
		h.SetMarkerColor( col + 1 )
		hists.append(h)

		#print group.name, h.Integral()/hists[0].Integral()

		if col== 0:
			l.AddEntry(h, group.name.split('_')[0])
		else:
			l.AddEntry(h, group.name.split('_')[0] + ' Reweighted')
		col += 1

	for i,h in enumerate(hists):
		if i == 0:
			h.SetMaximum(1.4*maxi)
			h.SetXTitle(entry.xtitle)
			h.Draw("PE")
		else:
			h.Draw("SAME,PE")

	l.Draw()



	tm = ROOT.TPaveText(0.19,0.74,0.22,0.88,"NDC")
	tm.SetBorderSize(0)
	tm.SetFillColor(0)
	tm.SetTextSize(0.048)
	tm.SetTextAlign(10)
	tmt = ''
	if 'mu' in mu.name:
	    tmt += '#mu#tau_{had}'
	elif 'el' in mu.name:
	    tmt += 'e#tau_{had}'
	tmt += '  '

	if hand == 'L':
		tmt += 'Left-Handed'
	elif hand == 'R':
		tmt += 'Right-Handed'

	tm.AddText(tmt)
	tm.AddText('#eta Reweighting')
	tm.Draw()


	c.cd(2)
	#ROOT.gPad.SetTopMargin(0)
	alpgen = hists[0].Clone('ratio')
	pythia = hists[1].Clone('rat1')
	pythia.Divide(alpgen)
	# powheg = hists[2].Clone('rat1')
	# powheg.Divide(alpgen)
	# herwig = hists[3].Clone('rat1')
	# herwig.Divide(alpgen)
	
	pythia.SetMaximum(1.08)
	pythia.SetMinimum(0.92)

	pythia.GetYaxis().SetNdivisions(5)
	pythia.GetYaxis().SetLabelSize(0.15)

	pythia.Draw("PE")
	# powheg.Draw("PE,SAME")
	# herwig.Draw("PE,SAME")
	# for i in xrange(1, len(hists)):

	# 	h = hists[i].Clone('rat_%i' % i)
	# 	h.Divide(alpgen)
	# 	if i == 1:
	# 		h.Draw("PE")
	# 		print h
	# 	else:
	# 		h.Draw("PE,SAME")
	#		print h

	#alpgen.Divide(alpgen)
	#alpgen.Draw("PE,SAME")
	xmin = alpgen.GetXaxis().GetXmin()
	xmax = entry.xmax or alpgen.GetXaxis().GetXmax()
	cl = ROOT.TLine(xmin,1,xmax,1)
	cl.SetLineStyle(3)

	cl.Draw("SAME")

	c.Update()

	c.Print('../plots/note/%s.eps' % cname)
	c.SaveAs('../plots/note/%s.gif+1' % pname)


def compare_truth ( mu, entry, hand, norm = False):

	cname = '%s_%s' % ( entry, hand)
	pname = 'Truth_plots' 
	if norm:
		cname += '_norm'


	c = ROOT.TCanvas(cname,cname,1500,1000)
	c.UseCurrentStyle()
	c.Divide(1,2)
	c.cd(1).SetPad(0.0,0.2,1.0,1.0)
	c.cd(2).SetPad(0.0,0.0,1.0,0.2)
	c.cd(1)

	l = ROOT.TLegend(.68,.7,.85,.9)
	l.SetTextSize(0.035)
	l.SetBorderSize(0)
	l.SetFillColor(0)


	if 'mZ' in entry:
		if 'L' in hand:
			histname = 'mZ_left_unw'
		else:
			histname = 'mZ_right_unw'
		c.SetLogy()
	elif not 'ups' in entry:
		histname = '%s_%s' % (entry, hand.lower())

	else:
		histname = '%s_%s_5' % (entry, hand.lower())
	hists = []
	maxi = 0
	col =0 
	for group in mu.MC:
		if not hand in group.name:
			continue
		col += 1

		gth = None
		for ds in group:	
			h_ds = None
			ds.load()
			th = None
			for f in ds.root_files:
				h = f.Get('massHists/%s' % histname)
				#hbin = h.GetNbinsX()
				#print hbin

			 	if not h_ds:
			 		h_ds = h
			 	else:
			 		h_ds.Add(h)
			h_ds.Scale( ds.xsec * int_lumi / ds.totalEvents)
			h_ds.SetLineColor(col)
			h_ds.SetMarkerColor(col)
			h_ds.SetXTitle(entry)
			h_ds.SetMarkerSize(1)
			h_ds.SetLineWidth(1)
			#ds.close()

			if not gth:
				gth = h_ds
			else:
				gth.Add(h_ds)

		hists.append(gth)
		maxi = max(gth.GetMaximum(),maxi)
		l.AddEntry(gth,group.name.split('_')[0])


	for i,h in enumerate(hists):
		if i == 0:
			h.SetMaximum(1.45*maxi)
			#h.SetXTitle(entry.xtitle)
			h.Draw("PE")
		else:
			h.Draw("SAME,PE")

	l.Draw()

	tm = ROOT.TPaveText(0.19,0.82,0.22,0.88,"NDC")
	tm.SetBorderSize(0)
	tm.SetFillColor(0)
	tm.SetTextSize(0.048)
	tm.SetTextAlign(10)
	tmt = ''
	if hand == 'L':
		tmt += 'Left-Handed'
	elif hand == 'R':
		tmt += 'Right-Handed'
	tm.AddText(tmt)
	tm.Draw()


	c.cd(2)
	#ROOT.gPad.SetTopMargin(0)
	alpgen = hists[0].Clone('ratio')
	pythia = hists[1].Clone('rat1')
	pythia.Divide(alpgen)
	powheg = hists[2].Clone('rat1')
	powheg.Divide(alpgen)
	herwig = hists[3].Clone('rat1')
	herwig.Divide(alpgen)
	
	pythia.SetMaximum(1.15)
	pythia.SetMinimum(0.85)

	pythia.GetYaxis().SetNdivisions(3)
	pythia.GetYaxis().SetLabelSize(0.15)

	pythia.Draw("P")
	powheg.Draw("P,SAME")
	herwig.Draw("P,SAME")
	# for i in xrange(1, len(hists)):

	# 	h = hists[i].Clone('rat_%i' % i)
	# 	h.Divide(alpgen)
	# 	if i == 1:
	# 		h.Draw("PE")
	# 		print h
	# 	else:
	# 		h.Draw("PE,SAME")
	#		print h

	#alpgen.Divide(alpgen)
	#alpgen.Draw("PE,SAME")
	xmin = alpgen.GetXaxis().GetXmin()
	xmax = alpgen.GetXaxis().GetXmax()
	cl = ROOT.TLine(xmin,1,xmax,1)
	cl.SetLineStyle(3)

	cl.Draw("SAME")

	#c.Update()

	c.SaveAs('../plots/note/%s.eps' % cname)
	c.SaveAs('../plots/note/%s.gif+1' % pname)


def truth_save( entry, hand ):
	rf = ROOT.TFile('TruthWeights.root','UPDATE')
	cname = '%s_%s' % ( entry, hand)

	if 'mZ' in entry:
		if 'L' in hand:
			histname = 'mZ_left_unw'
		else:
			histname = 'mZ_right_unw'
		c.SetLogy()
	elif not 'ups' in entry:
		histname = '%s_%s' % (entry, hand.lower())

	else:
		histname = '%s_%s_5' % (entry, hand.lower())
	hists = []
	maxi = 0
	col =0 
	for group in mu.MC:
		if not hand in group.name:
			continue
		col += 1

		gth = None
		for ds in group:	
			h_ds = None
			ds.load()
			th = None
			for f in ds.root_files:
				h = f.Get('massHists/%s' % histname)
				#hbin = h.GetNbinsX()
				#print hbin

			 	if not h_ds:
			 		h_ds = h
			 	else:
			 		h_ds.Add(h)
			h_ds.Scale( ds.xsec * int_lumi / ds.totalEvents)
			h_ds.SetLineColor(col)
			h_ds.SetMarkerColor(col)
			h_ds.SetXTitle(entry)
			h_ds.SetMarkerSize(1)
			h_ds.SetLineWidth(1)
			#ds.close()

			if not gth:
				gth = h_ds
			else:
				gth.Add(h_ds)

		gth.SetName('%s_%s_%s' % (group.name, entry, hand))
		rf.cd()
		gth.Write()


def truth_weight( entry, hand):
	rf = ROOT.TFile('TruthWeights.root','UPDATE')


	if 'mZ' in entry:
		if 'L' in hand:
			histname = 'mZ_left_unw'
		else:
			histname = 'mZ_right_unw'
		c.SetLogy()
	elif not 'ups' in entry:
		histname = '%s_%s' % (entry, hand.lower())

	else:
		histname = '%s_%s_5' % (entry, hand.lower())

	for group in mu.MC:
		if not hand in group.name:
			continue
		if 'Alpgen' in group.name:
			continue
		nth = '%s_%s_%s_%s' % ('AlpgenPythia', hand , entry, hand)
		gth = '%s_%s_%s' % (group.name, entry, hand)

		gen = rf.Get( gth  ).Clone('alpgen_%s' % group.name )
		nom = rf.Get( nth ).Clone(group.name)

		gen.Rebin(5)
		nom.Rebin(5)

		gen.Divide(nom)
		gen.SetName('weight_%s' % gth)
		gen.Write()


def draw_weights( entry, hand):
	rf = ROOT.TFile('TruthWeights.root')

	cname = 'weights_Pythia8_%s' % (entry)

	c = ROOT.TCanvas(cname,cname,1500,1000)
	c.UseCurrentStyle()
	#l = ROOT.TLegend(.68,.7,.85,.9)
	l = ROOT.TLegend(.70,.78,.9,.9)
	l.SetTextSize(0.035)
	l.SetBorderSize(0)
	l.SetFillColor(0)
		

	hand = 'L'
	h_L = rf.Get( 'weight_Pythia8_%s_%s_%s' % (hand, entry, hand))
	hand = 'R'
	h_R = rf.Get( 'weight_Pythia8_%s_%s_%s' % (hand, entry, hand))

	h_L.SetLineColor(2)
	h_R.SetLineColor(4)
	h_L.SetMarkerColor(2)
	h_R.SetMarkerColor(4)
	#h.SetTitle( 'Pythia8/AlpgenPythia' )
	l.AddEntry( h_L, 'Left-Handed')
	l.AddEntry(h_R, 'Right-Handed')
	h_R.SetXTitle( 'p_{T}' )
	h_R.SetMaximum(1.4*h_R.GetMaximum())
	h_R.SetMinimum(0.6*h_R.GetMinimum())
	h_R.Draw("")
	h_L.Draw("SAME")
	l.Draw()

	c.SaveAs('../plots/note/%s.eps' % cname)


def compare_polarization( mu, entry, gen, zoom = True ):
    cname = "Pol_Alpgen_%s" % gen
    if zoom:
    	cname += '_zoom'
    c = ROOT.TCanvas()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)

    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    hists = []
    maxi = 0
    mini = 0
    axis = None
    if entry in ['Zpt']:
        axis = [ 10*i for i in xrange(20)]
    if entry in ['mZ']:
    	if not zoom:
    		axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
    		axis.extend( [80 + i for i in xrange(19)] )
    		axis.extend(  [100 + 2*i for i in xrange (4) ])
    		axis.extend([110,115,120,125])
    		axis.extend([130+i*20 for i in xrange(3)])
    	else:
    		axis = [86 + 0.5*i for i in xrange(20)] 
        m_axis = []
     	for x in axis:
     		m_axis.append( x * 1000 )
     	axis = m_axis

   	hists = []
   	maxi = 0; mini = 0
   	col =0 
   	## alpgen
   	for i,group in enumerate(mu.MC):
   		if group.name.split('_')[0] not in ['AlpgenPythia', gen]:
   			continue


   		hand = group.name[-1]
   		if hand == 'L':
	   		gth_left = None
	   		if 'mZ' in entry:
				histname = 'mZ_left_unw'
   		elif hand == 'R':
	   		gth_right = None
	   		if 'mZ' in entry:
		   		histname = 'mZ_right_unw'
   		for ds in group:	
   			h_ds = None
   			ds.load()
   			th = None
   			for f in ds.root_files:
   				h = f.Get('massHists/%s' % histname)
   				#hbin = h.GetNbinsX()
   				#print hbin

   			 	if not h_ds:
   			 		h_ds = h
   			 	else:
   			 		h_ds.Add(h)

   			h_ds.Scale( ds.xsec * int_lumi * ds.xpol/ ds.totalEvents)
   			#print ds.name, ds.totalEvents
   			# h_ds.SetLineColor(col)
   			# h_ds.SetMarkerColor(col)
   			# h_ds.SetXTitle(entry)
   			# h_ds.SetMarkerSize(1)
   			# h_ds.SetLineWidth(1)
   			# #ds.close()

   			if 'L' in hand:
	   			if not gth_left:
   					gth_left = h_ds
   				else:
   					gth_left.Add(h_ds)
   			elif 'R' in hand:
   				if not gth_right:
   					gth_right = h_ds
   				else:
   					gth_right.Add(h_ds)

   		if i%2 == 1:
   			h = gth_right.Clone()
   			hT = gth_right.Clone()
   			h.Add(gth_left, -1)
   			hT.Add(gth_left)


   			if axis:
   			    h = h.Rebin(len(axis) -1 , "rebin_h_%s" % group.name, array('d',axis))
   			    hT = hT.Rebin(len(axis) -1 , "rebin_hT_%s" % group.name, array('d',axis))

   			else:
   			    try:
   			        h.Rebin( 6 )
   			        hT.Rebin( 6 )
   			    except:
   			        h.Rebin( 5)
   			        hT.Rebin( 5 )

   			#print h.Integral(), hT.Integral()

   			h.Divide(hT)

   			h.SetMarkerSize(0.3)
   			h.SetLineWidth(1)
   			if 'Alpgen' not in group.name:
   				h.SetMarkerColor(4)
   				h.SetLineColor(4)

   			maxi = max(maxi,h.GetMaximum())
   			mini = min(mini, h.GetMinimum())

   			hists.append(h)
   			l.AddEntry( h, group.name.split('_')[0]   )

   			#print h.Integral()




    draw = ''
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        h.SetMinimum( 1.4 * mini )
        if zoom:
        	h.GetXaxis().SetNdivisions(5)
        #h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    #tm.Draw() 
    l.Draw()

    c.SaveAs( '../plots/note/%s.eps' % (cname))    



def compare_pol_SR( mu, entry, gen, zoom = False ):
    #cname = "Pol_SR_%s_%s" % ("Reweight", mu.name)
    cname = "Pol_SR_%s_%s" % (gen, mu.name)
    if zoom:
        cname += '_zoom'

    rf = ROOT.TFile('~/panalysis/output/prod4/theory/%s_regions.root' % mu.name)

    c = ROOT.TCanvas()

    l = ROOT.TLegend(0.6,0.77,0.92,0.92)

    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    hists = []
    maxi = 0
    mini = 0
    axis = None
    if entry in ['Zpt']:
        axis = [ 10*i for i in xrange(20)]
    if False:
        if not zoom:
            axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
            axis.extend( [80 + i for i in xrange(19)] )
            axis.extend(  [100 + 2*i for i in xrange (4) ])
            axis.extend([110,115,120,125])
            #axis.extend([130+i*20 for i in xrange(3)])
        else:
            axis = [86 + 0.5*i for i in xrange(20)] 
        m_axis = []
        for x in axis:
            m_axis.append( x * 1000 )
        axis = m_axis

    hists = []
    maxi = 0; mini = 0
    col =0 
    region = 'SR_OS'
    ## alpgen


    #for i,group in enumerate(mu.MC[0:2]+mu.Alpgen):
    for i,group in enumerate(mu.MC):

        if group.name.split('_')[0] not in ['AlpgenPythia', gen]:
             continue

        print group.name, group.xpol

        hand = group.name[-1]
        if hand == 'L':
            gth_left = rf.Get('h_%s_%s_%s' % (entry.name, group.name, region))
            gth_left.Scale( group.xpol  )

        elif hand == 'R':
            gth_right = rf.Get('h_%s_%s_%s' % (entry.name, group.name, region))
            gth_right.Scale( group.xpol )


        if i%2 == 1:
            h = gth_right.Clone()
            hT = gth_right.Clone()
            h.Add(gth_left, -1)
            hT.Add(gth_left)


            if axis:
                h = h.Rebin(len(axis) -1 , "rebin_h_%s" % group.name, array('d',axis))
                hT = hT.Rebin(len(axis) -1 , "rebin_hT_%s" % group.name, array('d',axis))

            else:
                try:
                    h.Rebin( 6 )
                    hT.Rebin( 6 )
                except:
                    h.Rebin( 5)
                    hT.Rebin( 5 )

            #print h.Integral(), hT.Integral()

            h.Divide(hT)
            print h
            h.SetMarkerSize(0.3)
            h.SetLineWidth(1)
            if 'Pythia8' in group.name:
                h.SetMarkerColor(4)
                h.SetLineColor(4)

            maxi = max(maxi,h.GetMaximum())
            mini = min(mini, h.GetMinimum())

            hists.append(h)
            l.AddEntry( h, group.name.split('_')[0]   )
            #if not 'Pythia8' in group.name:
            #    l.AddEntry( h, group.name.split('_')[0]   )
            #else:
            #    l.AddEntry( h, group.name.split('_')[0]  + ' Reweighted'  )

            #print h.Integral()

    draw = ''
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        h.SetMinimum( 1.4 * mini )
        if zoom:
            h.GetXaxis().SetNdivisions(5)
        #h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'
        #break
    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    #tm.Draw() 
    l.Draw()

    c.SaveAs( '../plots/note/%s.eps' % (cname))    




if __name__ == '__main__':

	entry = trueTau_upsilon
	mu = streams.mutheory
	trueTau_upsilon.rebin = 2
	truethad_phi.rebin = 2
	truethad_eta.rebin = 2
	truethad_pt.rebin = 2

	# entry = 'hist_trueta'
	entries = ['hist_trueta','hist_truphi','hist_trupt','hist_truviseta','hist_truvisphi','hist_truvispt','hist_ups','mZ']
	entries = ['hist_trueta']
	entries = ['hist_trupt']
	for entry in entries:
	    #truth_save( entry, 'L'  )
	    #truth_save( entry, 'R'  ) 

	    #truth_weight( entry, 'L'  )
	    #truth_weight( entry, 'R'  ) 


	    #compare_truth(mu, entry, 'L')
	    #compare_truth(mu, entry, 'R')
	    draw_weights( entry, 'L')
	    draw_weights( entry, 'R')    

	    #break


