

function adjusting_TDS%(delat_pahse)
phase_XTDS=getfield(doocsread(app.add_phase), 'data');

x1=getfield(doocsread(app.add_COFM), 'data');
tem=doocswrite(app.add_phase, phase_XTDS+0.5);
x2=getfield(doocsread(app.add_COFM), 'data');
if x2-x1 >0
    phase_feed=3 ;% or -1
elseif x2-x1 <0
    phase_feed=-3;
end

%%%%%%%%
addr_TDS_amp='';
tem=doocswrite(app.add_phase, phase_XTDS)
x1=getfield(doocsread(app.add_COFM), 'data')
A0_TDS=getfield(doocsread(addr_TDS_amp), 'data');
x0=getfield(doocsread(app.add_COFM), 'data')
for A1=A0_TDS:5:71
    tem=doocswrite(addr_TDS_amp, A1);
    x2=getfield(doocsread(app.add_COFM), 'data');
    if x2>x1 
    while x2>x1
    tem=doocswrite(app.add_phase, phase_feed);
    x2=getfield(doocsread(app.add_COFM), 'data');
    if x2>2000
        phase_feed=-phase_feed
    end
    end
    end
    
    if x2<x1 
    while x2<x1
    tem=doocswrite(app.add_phase, phase_feed);
    x2=getfield(doocsread(app.add_COFM), 'data');
    if x2<50
        phase_feed=-phase_feed
    end
    end
    end
    
    
    

end
