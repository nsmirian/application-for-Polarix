function adjustingTDSFAST()
clear all
clc
%function adjusting_TDS%(delat_pahse)

ddr_cam        = 'FLASH.DIAG/CAMERA/OTR9FL2XTDS/';
    addr_charge='FLASH.DIAG/TOROID/7FL2XTDS/CHARGE.FLASH2';
    addr_TDSphase='FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
    add_COFM0='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.X.CM';
    add_COFM='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.X.MEAN';
    add2='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.X.MEAN';
    add_phase='FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
    addr_dipol='FLASH.MAGNETS/MAGNET.ML/D4FL2XTDS/CURRENT.SP';
    addr_amp='FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.AMPL';
    phase_XTDS=getfield(doocsread(add_phase), 'data');
x1=getfield(doocsread(add_COFM), 'data')
tem=doocswrite(add_phase, phase_XTDS+2);
pause(0.2)
x2=getfield(doocsread(add_COFM), 'data')

if x2-x1 >0
    phase_feed=0.5 ;% or -1
elseif x2-x1 <0
    phase_feed=-0.5;
end

%%%%%%%%
addr_TDS_amp=addr_amp;
tem=doocswrite(add_phase, phase_XTDS)
x1=getfield(doocsread(add_COFM), 'data')
A0_TDS=getfield(doocsread(addr_TDS_amp), 'data');
x0=getfield(doocsread(add_COFM), 'data')
for A1=A0_TDS:2:72
    tem=doocswrite(addr_TDS_amp, A1);
    pause(0.2)
    x2=getfield(doocsread(add_COFM), 'data');
 %   if x2>x0
%     while x2>x0
%      phase=getfield(doocsread(addr_TDSphase), 'data')
%     tem=doocswrite(add_phase, phase+phase_feed);
%     pause(0.5)
%     x2=getfield(doocsread(add_COFM), 'data');
    if x2>11
        phase_feed=-phase_feed
        while x2>11
        
        phase=getfield(doocsread(addr_TDSphase), 'data')
    tem=doocswrite(add_phase, phase+phase_feed);
    pause(0.5)
    x2=getfield(doocsread(add_COFM), 'data');
        end
    end
  %  end
   % end
    
    if x2<3
    while x2<3
        phase=getfield(doocsread(addr_TDSphase), 'data')
    tem=doocswrite(add_phase, phase_feed+phase);
    pause(0.5)
    x2=getfield(doocsread(add_COFM), 'data');
    %if x2<2
     %   phase_feed=-phase_feed
    end
    end
   % end
    
    
    

end
