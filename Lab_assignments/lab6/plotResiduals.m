function plotResiduals(est)
global dat

ref = [dat.f300, dat.f3000, dat.f30000, dat.ctrl];
name = ["300 pN/s","3000 pN/s","30000 pN/s","Neg control"];

figure;
for i = 1:size(ref,2)
    plot(abs(ref(:,i)-est(:,i)))
    hold on
end

xlabel('force rupture')
ylabel('absolute residual')
title('absolute residuals per condition')
legend(name)

end