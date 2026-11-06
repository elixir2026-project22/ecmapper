"""Worked example: hand-chosen layout for the iML1515 central-metabolism draft map (BRENDA map draft).

Input : examples/central_metabolism_escher.json   (91 reactions, 10 subsystems, node ids are those of that file)
Result: glycolysis as the central column; PDH + TCA ring below it; PPP and ED on the left; PEP/OAA node left of the
        spine; fermentations on the right; oxidative phosphorylation pieces packed into free space.

Read this file as a style reference for new specs: every number below is a hand-chosen coordinate.
Conventions: x grows to the right, y grows downwards, units are Escher map units (a primary node has radius 20).
"""
import numpy as np
from escher_handlayout.engine import unit


def build(L):
    # ------------------------------------------------------------------ subsystems: colour, title, existing title label
    L.zone('Glycolysis', '#2b5c9e', 'Glycolysis', 'glycolysis')
    L.zone('Pentose phosphate pathway', '#5b8c2a', 'Pentose phosphate pathway', 'pentose phosphate pathway')
    L.zone('Entner-Doudoroff pathway', '#8e3b9e', 'Entner-Doudoroff pathway', 'Entner Doudoroff pathway')
    L.zone('Lactate fermentation', '#d9822b', 'Lactate fermentation', 'lactate fermentation')
    L.zone('Oxidative decarboxylation of pyruvate', '#6a4fb3', 'Pyruvate dehydrogenase', 'oxidative decarboxylation of pyruvate')
    L.zone('Ethanol fermentation', '#8a7a2a', 'Ethanol fermentation / PFL, POR', 'ethanol fermentation')
    L.zone('Acetate fermentation', '#2f9a57', 'Acetate fermentation', 'acetate fermentation')
    L.zone('Citric acid cycle', '#c9a100', 'Citric acid cycle + glyoxylate shunt', 'citric acid cycle')
    L.zone('Gluconeogenesis', '#c2415a', 'PEP / pyruvate / OAA node', 'gluconeogenesis')
    L.zone('Oxidative phosphorylation', '#1a8f9a', 'Oxidative phosphorylation', 'oxidative phosphorylation')
    # titles for pieces that end up disconnected (used by the packing step): bigg_id of one of their reactions -> (zone, text, label to reuse)
    L.small_titles.update({
        'ATPM': ('Oxidative phosphorylation', 'oxidative phosphorylation: ATP synthase / ATP turnover', 'oxidative phosphorylation'),
        'NADH10': ('Oxidative phosphorylation', 'NADH dehydrogenase (menaquinone)', None),
        'NADH16pp': ('Oxidative phosphorylation', 'NADH dehydrogenase (ubiquinone)', None),
        'NADH18pp': ('Oxidative phosphorylation', 'NADH dehydrogenase (demethylmenaquinone)', None),
        'FE3DCITR5': ('Citric acid cycle', 'Fe(III)-citrate reduction', None),
        'FE3DCITexs': ('Citric acid cycle', 'Fe(III)-citrate formation (extracellular)', None),
        'PGM': ('Entner-Doudoroff pathway', 'ED: glycerate branch', None)})
    L.bend = 45
    # ------------------------------------------------------------------ node splits / merges
    xu_copy = L.split('871', 909)            # xylulose-5P for TKT2
    f6p_t2  = L.split('898', 909)            # TKT2 product F6P  (-> glycolysis F6P)
    accoa_m = L.split('147', 250)            # acetyl-CoA for malate synthase (inside TCA ring)
    cit_fe  = L.split('135', 981)            # citrate produced by Fe-citrate reductase (local copy)
    L.merge('480', '923')                    # glucose
    L.merge('482', '823')                    # G6P  glycolysis <-> PPP
    L.merge('493', f6p_t2)                   # F6P  glycolysis <-> PPP (TKT2)
    L.merge('528', '899')                    # G3P  glycolysis <-> PPP (TKT2)
    L.merge('604', '616', '656', '342', '414')   # pyruvate hub
    L.merge('592', '387')                    # PEP
    L.merge('657', '343', '147')             # acetyl-CoA hub (PDH/PFL/POR5/ACALD/CS)
    L.merge('355', '60')                     # acetaldehyde (ethanol <-> acetate fermentation)
    L.merge('149', '384')                    # oxaloacetate  TCA <-> PEP/OAA node
    L.merge('807', '1')                      # 6-phosphogluconate PPP <-> ED
    gx = 3000
    def pshift(a, b, off):
        pa, pb = np.array(L.pos[str(a)]), np.array(L.pos[str(b)]); u = unit(pb-pa); n = np.array([-u[1], u[0]])
        return tuple(n*off)

    # ================================================================== GLYCOLYSIS (central column)
    L.setp({480:(gx,500), 482:(gx,800), 493:(gx,1100), 504:(gx,1400), 527:(gx+300,1750), 528:(gx,1750),
            548:(gx,2100), 561:(gx,2400), 574:(gx,2700), 592:(gx,3000), 604:(gx,3750)})
    z='Glycolysis'
    L.reaction(492,z,side='R',dist=110)
    L.reaction(934,z,shift=(260,0),side='R',dist=110)
    L.reaction(1017,z,shift=(-110,0),side='L')
    L.reaction(501,z,side='R')
    L.reaction(514,z,side='R',dist=110)
    L.reaction(526,z,shift=(-110,0),side='L')
    L.reaction(537,z,axis=(0,1),mid=(gx,1580),lab='L')
    L.reaction(545,z,lab='D')
    L.reaction(560,z,side='R',dist=110)
    L.reaction(573,z,side='R',dist=110)
    L.reaction(582,z,shift=(-110,0),side='L')
    L.reaction(590,z,shift=(110,0),side='R')
    L.reaction(601,z,side='R',dist=110)
    for rid,off in ((1070,-225),(614,-75),(1082,75),(1094,225)):
        L.reaction(rid,z,shift=(0,off),side='R',curve=False,gap=40,dist=120)
    L.titles[z]=(gx+250,330)
    L.bend = 70

    # ================================================================== PENTOSE PHOSPHATE PATHWAY (left of glycolysis)
    z='Pentose phosphate pathway'
    L.setp({811:(400,800), 861:(900,1100), 871:(900,1750), 910:(1600,1100), 911:(1600,1750),
            897:(2300,1100), 898:(2300,1750), xu_copy:(2300,2400), 825:(2400,800), 807:(1650,800),
            870:(900,2300), 888:(400,2300), 1032:(1600,2350)})
    L.reaction(836,z,side='D')
    L.reaction(860,z,side='U')
    L.reaction(822,z,side='U')
    L.reaction(850,z,shift=(0,70),side='D')
    L.reaction(869,z,lab='R',bow=0.05)
    L.reaction(887,z,lab='L',bow=0.06)
    L.reaction(921,z,lab='R')
    L.reaction(944,z,lab='R')
    L.reaction(909,z,axis=(1,0),mid=(2650,1425),lab='D')
    L.reaction(1042,z,side='R')
    L.reaction(896,z,lab='U')
    L.reaction(879,z,lab='R')
    L.titles[z]=(-250,560)

    # ================================================================== ENTNER-DOUDOROFF (above PPP, branching at 6PGC)
    z='Entner-Doudoroff pathway'
    L.setp({2:(1650,300), 26:(1500,-200), 27:(1800,-200), 13:(2250,300), 48:(100,-200), 37:(700,-200), 38:(700,-700)})
    L.reaction(12,z,axis=(0,-1),side='L')
    L.reaction(36,z,axis=(0,-1),mid=(1650,150),lab='R')
    L.reaction(25,z,side='D')
    L.reaction(59,z,side='U')
    L.reaction(46,z,lab='R')
    L.titles[z]=(1000,-330)

    # ================================================================== LACTATE FERMENTATION (right of pyruvate)
    z='Lactate fermentation'
    L.setp({615:(4400,3250), 641:(4400,3650)})
    L.reaction(628,z,shift=pshift(615,604,-80),side='U',lab='D',dist=190)
    L.reaction(640,z,shift=pshift(615,604,80),side='D',lab='U',dist=190)
    L.reaction(653,z,side='D')
    L.titles[z]=(4350,2950)

    # ================================================================== PYRUVATE -> ACETYL-CoA (PDH, PFL, POR5)
    L.setp({657:(gx,4450)})
    z='Oxidative decarboxylation of pyruvate'
    L.reaction(671,z,side='R',dist=130)
    z2='Ethanol fermentation'
    L.reaction(354,z2,shift=(300,0),side='R',dist=120)
    L.reaction(1058,z2,shift=(620,0),side='R',dist=120)
    L.titles[z]=(3150,4850)

    # ================================================================== ETHANOL FERMENTATION (acetyl-CoA -> acetaldehyde -> ethanol)
    L.setp({355:(4700,4450), 370:(4900,3950), 1117:(5700,3900)})
    L.reaction(369,z2,side='D')
    L.reaction(382,z2,side='R')
    L.reaction(1127,z2,side='U',bow=0.05)
    L.titles[z2]=(5000,3600)

    # ================================================================== ACETATE FERMENTATION
    z='Acetate fermentation'
    dx,dy=2300,3750
    L.setp({63:(2400+dx,1300+dy), 119:(3050+dx,850+dy), 78:(2400+dx,1900+dy), 94:(1500+dx,1600+dy),
            946:(3200+dx,1300+dy), 945:(3200+dx,1900+dy), 994:(4000+dx,1900+dy), 995:(4000+dx,1300+dy),
            1096:(4000+dx,700+dy), 1097:(4700+dx,1300+dy)})
    L.reaction(75,z,side='L')
    L.reaction(134,z,side='U',dist=140,step=85)
    L.reaction(105,z,side='U')
    L.reaction(117,z,side='D')
    L.reaction(92,z,shift=(-120,0),side='L',step=70)
    L.reaction(956,z,lab='R')
    L.reaction(1005,z,lab='R')
    L.reaction(1107,z,axis=(0,1),mid=(4000+dx,1000+dy),side='L')
    L.titles[z]=(5100,6000)

    # ================================================================== TCA CYCLE (below the acetyl-CoA hub)
    z='Citric acid cycle'
    C=np.array((3000,6600)); R_=1050
    ang = {135:-90, 136:-50, 163:-10, 223:30, 265:70, 251:110, 278:150, 173:190, 149:230}
    for n,a in ang.items():
        L.pos[str(n)] = (C[0]+R_*np.cos(np.radians(a)), C[1]+R_*np.sin(np.radians(a)))
    L.setp({accoa_m:(2750,6050), 237:(3100,6900), 1108:(4250,5850), 957:(2300,5150),
            967:(300,7300), 969:(1100,7550), cit_fe:(1100,7050), 982:(-200,7950), 983:(600,7950), 984:(200,7650)})
    def lane(a,b,off):
        pa,pb=np.array(L.pos[str(a)]),np.array(L.pos[str(b)]); m=(pa+pb)/2; n=(m-C)/np.linalg.norm(m-C)
        return tuple(n*off)
    cen=tuple(C)
    L.reaction(161,z,arc=(cen,0),arc_nodes=('149','135'),side='out')
    L.reaction(146,z,arc=(cen,0),side='out')
    L.reaction(172,z,arc=(cen,0),side='out')
    L.reaction(236,z,arc=(cen,0),side='out')
    L.reaction(277,z,arc=(cen,0),side='out')
    L.reaction(304,z,arc=(cen,0),side='out')
    L.reaction(288,z,arc=(cen,0),side='out')
    L.reaction(328,z,arc=(cen,-150),side='in')
    L.reaction(340,z,arc=(cen,0),side='out')
    L.reaction(316,z,arc=(cen,150),side='out')
    L.reaction(221,z,arc=(cen,-300),side='in')
    L.reaction(209,z,arc=(cen,-150),side='in')
    L.reaction(185,z,arc=(cen,0),side='out')
    L.reaction(197,z,arc=(cen,150),side='out')
    L.reaction(1029,z,arc=(cen,300),side='out')
    pi_=np.array(L.pos['163']); cp=(np.array(L.pos['237'])+np.array(L.pos['251']))/2; vec=cp-pi_
    L.reaction(260,z,axis=tuple(vec),mid=tuple(pi_+0.27*vec),side='U',bow=-0.08)
    L.reaction(250,z,side='D',lab='D')
    L.reaction(966,z,side='U')
    L.reaction(1116,z,lab='U')
    L.reaction(981,z,side='U')
    L.reaction(993,z,lab='R')
    L.titles[z]=(3500,7900)

    # ================================================================== PEP / OAA node (gluconeogenesis), left of the spine
    z='Gluconeogenesis'
    L.setp({454:(1900,3450), 439:(1300,4150)})
    L.reaction(398,z,shift=(-70,0),side='L',lab='L')
    L.reaction(411,z,shift=(70,0),side='R',lab='R')
    L.reaction(428,z,shift=(-170,0),side='L',step=70)
    L.reaction(438,z,side='R',bow=0.05)
    L.reaction(453,z,side='D',bow=0.05)
    L.reaction(465,z,shift=pshift(454,604,-70),side='U')
    L.reaction(478,z,shift=pshift(454,604,70),side='D')
    L.titles[z]=(700,3000)

    # ================================================================== OXIDATIVE PHOSPHORYLATION (bottom right)
    z='Oxidative phosphorylation'
    ox=6200; oy=6900
    L.setp({672:(ox-600,oy), 674:(ox+600,oy), 686:(ox-600,oy+470), 687:(ox+600,oy+470),
            708:(ox-600,oy+940), 709:(ox+600,oy+940)})
    L.reaction(685,z,shift=(0,-45),side='U'); L.reaction(707,z,shift=(0,45),side='D')
    L.reaction(697,z,shift=(0,-45),side='U'); L.reaction(729,z,shift=(0,45),side='D')
    L.reaction(719,z,shift=(0,-45),side='U'); L.reaction(739,z,shift=(0,45),side='D')
    hx=8300
    L.setp({740:(hx,7300), 775:(hx,7950)})
    L.reaction(786,z,shift=(-190,0),side='L')
    L.reaction(796,z,side='R',gap=85,secdist=430,curve=True,lab='R')
    L.reaction(806,z,shift=(190,0),side='R')
    L.reaction(753,z,axis=(0,-1),mid=(hx,6900),shift=(-190,0),side='L')
    L.reaction(763,z,axis=(0,1),mid=(hx,6900),side='R',gap=85,secdist=430,lab='R')
    L.reaction(773,z,axis=(0,-1),mid=(hx,6900),shift=(190,0),side='R')
    L.titles[z]=(5450,6500)
