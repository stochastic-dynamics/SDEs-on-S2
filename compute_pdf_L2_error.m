% compute_pdf_L2_error.m
%
% Compute the L2 error between two 2‑D PDFs and report it *as a percentage*
% of the L2 norm of the reference PDF (global relative error).  The script
% also plots the point‑wise squared error surface.
%
% -------- Usage --------
%   compute_pdf_L2_error();                        %% default file names
%   compute_pdf_L2_error('ref.mat','test.mat');    %% custom files
% -------------------------------------------------------------------------
function compute_pdf_L2_error( benchSamples, testSamples )

% ── 0.  Input arguments & defaults ──────────────────────────────────────
% if nargin < 1 || isempty(refFile)
%     refFile  = 'Basin_sig1_benchmark.mat';
% end
% if nargin < 2 || isempty(testFile)
%     testFile = 'Basin_sig1_weak1.mat';
% end
% if ~isfile(refFile) || ~isfile(testFile)
%     error('Cannot find %s or %s in the current folder.',refFile,testFile);
% end
% 
% % ── 1.  Load samples ────────────────────────────────────────
% ref  = load(refFile );
% test = load(testFile);
% 
% % 1.1  Extract [q, qdot] pairs  (edit to match your variable names)
% if isfield(ref,'xx') && isfield(test,'xx')
%     benchSamples = ref.xx;
%     testSamples  = test.xx;
% elseif isfield(ref,'qval1') && isfield(ref,'qdval1')
%     benchSamples = [ref.qval1(:)  ref.qdval1(:) ];
%     testSamples  = [test.qval1(:) test.qdval1(:)];
% else
%     error(['Sample arrays not found.  Please edit Section 1.1 of ',mfilename,'.m']);
% end
%
benchSamples = benchSamples(all(isfinite(benchSamples),2),:);
testSamples  = testSamples (all(isfinite(testSamples ),2),:);
%
testSamples(testSamples<-100) = 0;
testSamples(testSamples>100) = 0;
% ── 2.  Common evaluation grid ──────────────────────────────────────────
ngrid = 200;                            % resolution per axis
qGrid  = linspace(-1, 1, ngrid);
qdGrid = linspace(-5, 5, ngrid);
[Q,QD] = meshgrid(qGrid,qdGrid);
pts    = [Q(:) QD(:)];                 % N×2 query points
%
% ── 3.  Kernel‑density estimates on that grid ───────────────────────────
[~,~,bw] = ksdensity(benchSamples,'Support','unbounded');
if isscalar(bw); bw = [bw bw]; end     % ensure 1×2 bandwidth vector
%
fRef  = ksdensity(benchSamples, pts, 'Bandwidth',bw,'Support','unbounded');
fTest = ksdensity(testSamples , pts, 'Bandwidth',bw,'Support','unbounded');
%
fRef  = reshape(fRef ,ngrid,ngrid);
fTest = reshape(fTest,ngrid,ngrid);
%
dx = qGrid(2)-qGrid(1);
dy = qdGrid(2)-qdGrid(1);
%
% ── 4.  Global errors ─────────────────────────────────────────────────--
L2_abs = sqrt( sum( (fRef(:)-fTest(:)).^2 ) * dx * dy );   % absolute
L2_ref = sqrt( sum(  fRef(:)       .^2 ) * dx * dy );               % norm of ref
perc_global_err = 100 * L2_abs / L2_ref;                           % percentage
%
fprintf('\nL2 absolute error         : %.6g\n', L2_abs);
fprintf('L2 norm of PDF_{ref}     : %.6g\n', L2_ref);
fprintf('Percentage global error  : %.4f %%\n', perc_global_err);
%
% ── 5.  Plot squared error surface ─────────────────────────────────────
figure('Color',[1 1 1]);
imagesc(qGrid,qdGrid,(fRef-fTest).^2);
axis xy tight; set(gca,'YDir','normal'); colorbar;
xlabel('q','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('$\dot{\bf{q}}$','FontWeight','bold','FontName','Times New Roman','Interpreter','latex');
%
% if verLessThan('matlab','9.13')
    title('Point-wise squared error');
% else
    % title('Point-wise squared error');
% end
grid on; drawnow;
%
%
end  % function compute_pdf_L2_error
%
% End