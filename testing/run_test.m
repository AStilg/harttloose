% ultimate_spheroidal_formulation_test - is a test suite for several of the basic functions.
%
% NOTES: Originally to determine if the wavefunctions are appropriately locally orthogonal.
% It has some chance of failing this test. However, re-running the code should make it complete.
% If it does not, there are problems.
addpath ../ ../provided/ ./needed/

warn_state=false;
%% test 1/1a
batchno=1;
testnum=1;
disp(['------------BATCH ' num2str(batchno) '-----------']);
batchno=batchno+1;
% BATCH 1 tests important coordinate conversions. if these don't match then
% we can't be sure that the coordinate system is valid.

% verify that the conversion between coordinates is accurate:
J=[0.4337208491287741	-0.1327993463965178	-0.8912073600614353
0.8521572268758445	-0.2609187983088673	0.4535961214255774
0.2927700218845600	0.9561828874675149	0]';

cpo=1.3;
xepp=[2,1/3,11/10];
[Jtest,xyzp]=xietaphiv2xyzv(true,cpo,eye(3),[xepp;xepp;xepp]);

if sqrt(sum((J-Jtest).^2,'all'))>1e-14
    warning(['test ' num2str(testnum) 'a: prolate xietaphiv2xyzv is not working'])
else
    disp(['test ' num2str(testnum) 'a: prolate xietaphiv2xyzv passsed.'])
end
xzyr=[0.9629351738917601,1.891936182208246,0.8666666666666667];
if sqrt(sum((xyzp(1,:)-xzyr).^2))>1e-14
    warning(['test ' num2str(testnum) 'b: prolate xietaphi2xyz is not working'])
else
    disp(['test ' num2str(testnum) 'b: prolate xietaphi2xyz passsed.'])
end

[Jeye,xep2]=xyzv2xietaphiv(true,cpo,Jtest,xyzp);

if sqrt(sum((Jeye-eye(3)).^2,'all'))>1e-14
    warning(['test ' num2str(testnum) 'c: prolate xyzv2xietaphiv is not working'])
else
    disp(['test ' num2str(testnum) 'c: prolate xyzv2xietaphiv passsed.'])
end

if sqrt(sum((xep2-xepp).^2))>1e-14
    warning(['test ' num2str(testnum) 'd: prolate xyz2xietaphi is not working'])
else
    disp(['test ' num2str(testnum) 'd: prolate xyz2xietaphi passsed.'])
end

J=[0.4057087047351359	-0.2028543523675680	-0.8912073600614353
0.7971200956582002	-0.3985600478291001	0.4535961214255774
0.4472135954999579	0.8944271909999159	0]';

xepo=[1,1/3,11/10];
[Jtest,xyzo]=xietaphiv2xyzv(false,cpo,eye(3),[xepo;xepo;xepo]);

if sqrt(sum((J-Jtest).^2,'all'))>1e-14
    warning(['test ' num2str(testnum) 'e: oblate xietaphiv2xyzv is not working'])
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'e: oblate xietaphiv2xyzv passsed.'])
end

xzyr=[0.7862332771376675,1.544759424106488,0.4333333333333333];
if sqrt(sum((xyzo(1,:)-xzyr).^2))>1e-14
    warning(['test ' num2str(testnum) 'f: oblate xietaphi2xyz is not working'])
else
    disp(['test ' num2str(testnum) 'f: oblate xietaphi2xyz passsed'])
end

[Jeye,xep2]=xyzv2xietaphiv(false,cpo,Jtest,xyzo);

if sqrt(sum((Jeye-eye(3)).^2,'all'))>1e-14
    warning(['test ' num2str(testnum) 'g: oblate xyzv2xietaphiv is not working'])
else
    disp(['test ' num2str(testnum) 'g: oblate xyzv2xietaphiv passsed.'])
end

if sqrt(sum((xep2-xepo).^2))>1e-14
    warning(['test ' num2str(testnum) 'h: oblate xyz2xietaphi is not working'])
else
    disp(['test ' num2str(testnum) 'h: oblate xyz2xietaphi passsed.'])
end
testnum=testnum+1;
%% prolate test 2
disp(['------------BATCH ' num2str(batchno) '-----------'])
batchno=batchno+1;
eta=linspace(-1,1,9);
n0=0+[0:2:2*10];
m0=0;
S1ref1=[[0.2180710664879469,0.2574274193534938,0.2879909556157383,0.3073597286424433,0.3139929090466117,0.3073597286424433,0.2879909556157383,0.2574274193534938,0.2180710664879469];[0.6218207580992870,0.2565315027931586,-0.03904703122658538,-0.2310949605313757,-0.2976477390275101,-0.2310949605313757,-0.03904703122658538,0.2565315027931586,0.6218207580992870];[0.8457382317135784,-0.2791950668398425,-0.2557049323410043,0.1224577011305248,0.3126621380740191,0.1224577011305248,-0.2557049323410043,-0.2791950668398425,0.8457382317135784];[1.016971813721233,-0.2972590483406774,0.3238396256737235,0.03140195711089613,-0.3156672816921719,0.03140195711089613,0.3238396256737235,-0.2972590483406774,1.016971813721233];[1.163055441710774,0.2192484610470875,-0.07583919943145433,-0.1810538562903630,0.3167777014635734,-0.1810538562903630,-0.07583919943145433,0.2192484610470875,1.163055441710774];[1.292696688261680,0.3468417557498750,-0.2485579320001120,0.2872481755565094,-0.3173098694841252,0.2872481755565094,-0.2485579320001120,0.3468417557498750,1.292696688261680];[1.410461007231105,-0.1381592864992080,0.3275654857079593,-0.3226863092144080,0.3176058603712512,-0.3226863092144080,0.3275654857079593,-0.1381592864992080,1.410461007231105];[1.519119280996756,-0.3799381396181944,-0.08120283409961459,0.2780468979557040,-0.3177874005719945,0.2780468979557040,-0.08120283409961459,-0.3799381396181944,1.519119280996756];[1.620506344431256,0.04629046153674626,-0.2462545803081791,-0.1641613408898968,0.3179067520579985,-0.1641613408898968,-0.2462545803081791,0.04629046153674626,1.620506344431256];[1.715912381276767,0.3912602923635486,0.3286132284459187,0.009268168202905450,-0.3179894053146979,0.009268168202905450,0.3286132284459187,0.3912602923635486,1.715912381276767];[1.806285734556458,0.04958362951973630,-0.08337629235671851,0.1480588413313712,0.3180490055204672,0.1480588413313712,-0.08337629235671851,0.04958362951973630,1.806285734556458]];
isProlate=true;
S1test1=spheroidalS1(isProlate,n0(end),m0,1.5,eta,true);

if any(sqrt(sum(abs(S1ref1-(-1).^m0*S1test1(1:size(S1ref1,1),:)).^2))>1e-13)
    warning(['test ' num2str(testnum) 'a: prolate spheroidalS1 failed test, n0=0, m0=0 case 1.'])
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'a: prolate spheroidalS1 passed test case 1.'])
end
   
eta=linspace(-1,1,9);
n0=2+[0:2:2*10];
m0=1;
S1ref2=[[0,-0.3658200682673630,-0.3518903955279606,-0.2083181116050839,0,0.2083181116050839,0.3518903955279606,0.3658200682673630,0];[0,-0.2606928343599532,0.2269066719523473,0.2949543785858294,0,-0.2949543785858294,-0.2269066719523473,0.2606928343599532,0];[0,0.2701429756516749,0.1019637214344430,-0.3198716531841612,0,0.3198716531841612,-0.1019637214344430,-0.2701429756516749,0];[0,0.3183674376415376,-0.3328873166757278,0.2664771195929801,0,-0.2664771195929801,0.3328873166757278,-0.3183674376415376,0];[0,-0.1912966849871553,0.2356744274834902,-0.1467510912690510,0,0.1467510912690510,-0.2356744274834902,0.1912966849871553,0];[0,-0.3639239389572474,0.09543407196577798,-0.009872287433288156,0,0.009872287433288156,-0.09543407196577798,0.3639239389572474,0];[0,0.1012073426313557,-0.3318994342807897,0.1644401662311755,0,-0.1644401662311755,0.3318994342807897,-0.1012073426313557,0];[0,0.3884479764031849,0.2379727141672284,-0.2783556626092246,0,0.2783556626092246,-0.2379727141672284,-0.3884479764031849,0];[0,-0.004766882445550729,0.09318558564004648,0.3230894702721497,0,-0.3230894702721497,-0.09318558564004648,0.004766882445550729,0];[0,-0.3893733317334141,-0.3314740762959146,-0.2873541719450335,0,0.2873541719450335,0.3314740762959146,0.3893733317334141,0];[0,-0.09209406436963198,0.2390249305865793,0.1799591129035048,0,-0.1799591129035048,-0.2390249305865793,0.09209406436963198,0]];
isProlate=true;
S1test2=spheroidalS1(isProlate,n0(end),m0,2.1,eta,true);

if any(sqrt(sum(abs(S1ref2-(-1).^m0*S1test2(1:size(S1ref2,1),:)).^2))>1e-13)
    warning(['test ' num2str(testnum) 'b: prolate spheroidalS1 failed test, n0=2, m0=1 case 2.'])
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'b: prolate spheroidalS1 passed test case 2.'])
end

eta=linspace(-1,1,9);
n0=3+[0:2:2*10];
m0=3;
S1ref3=[[0,0.1157555474366588,0.2674883434847847,0.3803831708651249,0.4214793658262915,0.3803831708651249,0.2674883434847847,0.1157555474366588,0];[0,0.4035704137975370,0.2883360876049856,-0.1319626294240135,-0.3450074130864848,-0.1319626294240135,0.2883360876049856,0.4035704137975370,0];[0,0.3626230738265020,-0.3242328029559619,-0.06096879626710031,0.3310459798630810,-0.06096879626710031,-0.3242328029559619,0.3626230738265020,0];[0,-0.1272996826240365,-0.003113162518863638,0.2185001752156593,-0.3258599734829028,0.2185001752156593,-0.003113162518863638,-0.1272996826240365,0];[0,-0.4058459145402250,0.3119875032698481,-0.3116416807040778,0.3233311273897477,-0.3116416807040778,0.3119875032698481,-0.4058459145402250,0];[0,-0.03885900967629822,-0.2849456942677765,0.3214050492010009,-0.3218992397608405,0.3214050492010009,-0.2849456942677765,-0.03885900967629822,0];[0,0.3855927266928728,-0.03772563968000258,-0.2479972084514910,0.3210069042285412,-0.2479972084514910,-0.03772563968000258,0.3855927266928728,0];[0,0.1660257343488087,0.3184570198746807,0.1116049428098196,-0.3204121272329067,0.1116049428098196,0.3184570198746807,0.1660257343488087,0];[0,-0.3309608593710987,-0.2712955201640509,0.05243732656948878,0.3199953521659224,0.05243732656948878,-0.2712955201640509,-0.3309608593710987,0];[0,-0.2646231298797569,-0.05220809076358776,-0.2023943507799513,-0.3196917640945698,-0.2023943507799513,-0.05220809076358776,-0.2646231298797569,0];[0,0.2526544177070401,0.3215406845460179,0.3004873057700821,0.3194636646822254,0.3004873057700821,0.3215406845460179,0.2526544177070401,0]];
isProlate=true;
S1test3=spheroidalS1(isProlate,n0(end),m0,1.3,eta,true);

if any(sqrt(sum(abs(S1ref3-(-1).^m0*S1test3(1:size(S1ref2,1),:)).^2))>1e-13)
    warning(['test ' num2str(testnum) 'c: prolate spheroidalS1 failed test, n0=3, m0=3 case 3.'])
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'c: prolate spheroidalS1 passed test case 3.'])
end
testnum=testnum+1;

%% oblate test 2
eta=linspace(-1,1,9);
n0=0+[0:2:2*10];
m0=0;
S1ref1=[[0.3576572261459682,0.3047122243542183,0.2699640188583711,0.2502822938640426,0.2439103740367759,0.2502822938640426,0.2699640188583711,0.3047122243542183,0.3576572261459682];[0.6276758206598160,0.1694886143551662,-0.1213653502395404,-0.2821395342185245,-0.3335215126870512,-0.2821395342185245,-0.1213653502395404,0.1694886143551662,0.6276758206598160];[0.8463908992490861,-0.3117121173993015,-0.2322722043145171,0.1447907230643906,0.3221939611986972,0.1447907230643906,-0.2322722043145171,-0.3117121173993015,0.8463908992490861];[1.017140278656056,-0.2733312495382539,0.3332239239286081,0.01783827532377595,-0.3200552003783659,0.01783827532377595,0.3332239239286081,-0.2733312495382539,1.017140278656056];[1.163119917742866,0.2401484703583295,-0.09544477138255508,-0.1734944822705901,0.3193063271224713,-0.1734944822705901,-0.09544477138255508,0.2401484703583295,1.163119917742866];[1.292727157084442,0.3364341678905838,-0.2379370919753480,0.2845907749330002,-0.3189553863611458,0.2845907749330002,-0.2379370919753480,0.3364341678905838,1.292727157084442];[1.410477467717203,-0.1544180494467578,0.3317153292215364,-0.3238532392308319,0.3187624269112245,-0.3238532392308319,0.3317153292215364,-0.1544180494467578,1.410477467717203];[1.519129039464368,-0.3759252207369342,-0.09249154969136422,0.2817649101879623,-0.3186448693074153,0.2817649101879623,-0.09249154969136422,-0.3759252207369342,1.519129039464368];[1.620512539131450,0.05944819772507788,-0.2393698226686901,-0.1690581316465683,0.3185679099588035,-0.1690581316465683,-0.2393698226686901,0.05944819772507788,1.620512539131450];[1.715916525590798,0.3910488060363537,0.3312584496507292,0.01405847260259903,-0.3185147664171221,0.01405847260259903,0.3312584496507292,0.3910488060363537,1.715916525590798];[1.806288624834624,0.03896123711051125,-0.09131355574682244,0.1443881110281040,0.3184765226569858,0.1443881110281040,-0.09131355574682244,0.03896123711051125,1.806288624834624]];
isProlate=false;
S1test1=spheroidalS1(isProlate,n0(end),m0,1.5,eta,true);

if any(sqrt(sum(abs(S1ref1-(-1).^m0*S1test1(1:size(S1ref1,1),:)).^2))>1e-13)
    warning(['test xa: oblate spheroidalS1 failed test, n0=0, m0=0 case 1.'])
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'a: oblate spheroidalS1 passed test case 1.'])
end

eta=linspace(-1,1,9);
n0=2+[0:2:2*10];
m0=1;
S1ref2=[[0.*10^-67,-0.3978602001327477,-0.3144209616685084,-0.1654088633163745,0.*10^-67,0.1654088633163745,0.3144209616685084,0.3978602001327477,0.*10^-67];[0.*10^-67,-0.1760646553308249,0.2833789598176307,0.2906068129992267,0.*10^-67,-0.2906068129992267,-0.2833789598176307,0.1760646553308249,0.*10^-67];[0.*10^-67,0.3134929244979045,0.05312999285186786,-0.3286599180292792,0.*10^-67,0.3286599180292792,-0.05312999285186786,-0.3134929244979045,0.*10^-67];[0.*10^-67,0.2851379728865079,-0.3243107763382616,0.2810981046215691,0.*10^-67,-0.2810981046215691,0.3243107763382616,-0.2851379728865079,0.*10^-67];[0.*10^-67,-0.2256917221473357,0.2585572235779588,-0.1626484246325579,0.*10^-67,0.1626484246325579,-0.2585572235779588,0.2256917221473357,0.*10^-67];[0.*10^-67,-0.3492849049002804,0.07025572770853582,0.003992131400727121,0.*10^-67,-0.003992131400727121,-0.07025572770853582,0.3492849049002804,0.*10^-67];[0.*10^-67,0.1294132098864273,-0.3265779929979871,0.1547710195778591,0.*10^-67,-0.1547710195778591,0.3265779929979871,-0.1294132098864273,0.*10^-67];[0.*10^-68,0.3838438719841369,0.2523116853668560,-0.2738543649390767,0.*10^-68,0.2738543649390767,-0.2523116853668560,-0.3838438719841369,0.*10^-68];[0.*10^-68,-0.02791673638909548,0.07618138212675985,0.3235962685327823,0.*10^-68,-0.3235962685327823,-0.07618138212675985,0.02791673638909548,0.*10^-68];[0.*10^-68,-0.3910312613440520,-0.3276110329517411,-0.2917851890783259,0.*10^-68,0.2917851890783259,0.3276110329517411,0.3910312613440520,0.*10^-68];[0.*10^-68,-0.07345377569047437,0.2494682089474369,0.1866416903812713,0.*10^-68,-0.1866416903812713,-0.2494682089474369,0.07345377569047437,0.*10^-68]];
isProlate=false;
S1test2=spheroidalS1(isProlate,n0(end),m0,2.1,eta,true);

if any(sqrt(sum(abs(S1ref2-(-1).^m0*S1test2(1:size(S1ref2,1),:)).^2))>1e-13)
    warning(['test ' num2str(testnum) 'b: oblate spheroidalS1 failed test, n0=2, m0=1 case 2.'])
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'b: oblate spheroidalS1 passed test case 2.'])
end

eta=linspace(-1,1,9);
n0=3+[0:2:2*10];
m0=3;
S1ref3=[[0.*10^-68,0.1259947243602795,0.2745547617250671,0.3769218409440665,0.4127696878652876,0.3769218409440665,0.2745547617250671,0.1259947243602795,0.*10^-68];[0.*10^-69,0.4096351461305913,0.2731249446543309,-0.1429504705020929,-0.3470172773376683,-0.1429504705020929,0.2731249446543309,0.4096351461305913,0.*10^-69];[0.*10^-69,0.3492770510876822,-0.3313446375320039,-0.05273960613096201,0.3327801797571982,-0.05273960613096201,-0.3313446375320039,0.3492770510876822,0.*10^-69];[0.*10^-69,-0.1425533117699792,0.009764342763543111,0.2139826555588123,-0.3271024509280661,0.2139826555588123,0.009764342763543111,-0.1425533117699792,0.*10^-69];[0.*10^-70,-0.4038735354531117,0.3073385955628783,-0.3105064913538628,0.3242344301108680,-0.3105064913538628,0.3073385955628783,-0.4038735354531117,0.*10^-70];[0.*10^-70,-0.02700167645806250,-0.2903556520256066,0.3228824136489824,-0.3225778648772161,0.3228824136489824,-0.2903556520256066,-0.02700167645806250,0.*10^-70];[0.*10^-70,0.3880569173850770,-0.02978032751045615,-0.2510888657985964,0.3215327480862454,-0.2510888657985964,-0.02978032751045615,0.3880569173850770,0.*10^-70];[0.*10^-70,0.1575255348104404,0.3158343173748936,0.1152619306546548,-0.3208304614550362,0.1152619306546548,0.3158343173748936,0.1575255348104404,0.*10^-70];[0.*10^-70,-0.3354021005847465,-0.2753222052961837,0.04914449395216046,0.3203355713536300,0.04914449395216046,-0.2753222052961837,-0.3354021005847465,0.*10^-70];[0.*10^-70,-0.2589059053475627,-0.04649201569653835,-0.2001385078180062,-0.3199736123530177,-0.2001385078180062,-0.04649201569653835,-0.2589059053475627,0.*10^-70];[0.*10^-70,0.2579078256228564,0.3197372042183245,0.2996042633778582,0.3197008294554396,0.2996042633778582,0.3197372042183245,0.2579078256228564,0.*10^-70]];
isProlate=false;
S1test3=spheroidalS1(isProlate,n0(end),m0,1.3,eta,true);

if any(sqrt(sum(abs(S1ref3-(-1).^m0*S1test3(1:size(S1ref2,1),:)).^2))>1e-13)
    warning(['test ' num2str(testnum) 'c: oblate spheroidalS1 failed test, n0=3, m0=3 case 3.'])
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'c: oblate spheroidalS1 passed test case 3.'])
end
testnum=testnum+1;
%%
disp(['------------BATCH ' num2str(batchno) '-----------'])
batchno=batchno+1;
%random selection of modes;
n_modes=10;
n=zeros(n_modes,1);
m=n;
n=randi(11,n_modes,1)-1;
for ii=1:length(n)
    m(ii)=randi(n(ii)+1)-1;
end
m=sign(randn(n_modes,1)).*m;
%% test below...

isProlate=true;
c=abs(randn);

x=randn+linspace(-1e-8,1e-8,3);
y=randn+linspace(-1e-8,1e-8,3);
z=randn+linspace(-1e-8,1e-8,3);
[X,Y,Z]=meshgrid(x,y,z);

xep=xyz2xietaphi(isProlate,c,[X(:),Y(:),Z(:)]);
[M,N]=spheroidalvwf(isProlate,n,m,c,xep(:,1),xep(:,2),xep(:,3),3);
%

Mfull=reshape(M.',[prod(size(X)),3,n_modes]);
Nfull=reshape(N.',[prod(size(X)),3,n_modes]);


Mx=zeros(3,3,3);
My=Mx;
Mz=My;
MX=zeros(3,3,3,n_modes);
MY=MX;
MZ=MX;
NX=MX;
NY=MX;
NZ=MX;
CMX=zeros(3,3,3,n_modes);
CMY=CMX;
CMZ=CMX;
CNX=CMX;
CNY=CMY;
CNZ=CMZ;

for jj=1:n_modes
    %extract into cart.
    [Mx(:),My(:),Mz(:)]=xietaphiv2xyzv(isProlate,c,Mfull(:,:,jj),xep(:,1:3));
    MX(:,:,:,jj)=Mx;
    MY(:,:,:,jj)=My;
    MZ(:,:,:,jj)=Mz;
    [CMX(:,:,:,jj),CMY(:,:,:,jj),CMZ(:,:,:,jj)]=curl(X,Y,Z,Mx,My,Mz);
    [Mx(:),My(:),Mz(:)]=xietaphiv2xyzv(isProlate,c,Nfull(:,:,jj),xep(:,1:3));
    NX(:,:,:,jj)=Mx;
    NY(:,:,:,jj)=My;
    NZ(:,:,:,jj)=Mz;
    [CNX(:,:,:,jj),CNY(:,:,:,jj),CNZ(:,:,:,jj)]=curl(X,Y,Z,Mx,My,Mz);
end

%%


isProlate=false;
c=abs(randn);

x=randn+linspace(-1e-8,1e-8,3);
y=randn+linspace(-1e-8,1e-8,3);
z=randn+linspace(-1e-8,1e-8,3);
[X,Y,Z]=meshgrid(x,y,z);

xep=xyz2xietaphi(isProlate,c,[X(:),Y(:),Z(:)]);
[M,N]=spheroidalvwf(isProlate,n,m,c,xep(:,1),xep(:,2),xep(:,3),3);
%

Mfull=reshape(M.',[prod(size(X)),3,n_modes]);
Nfull=reshape(N.',[prod(size(X)),3,n_modes]);


Mx=zeros(3,3,3);
My=Mx;
Mz=My;
MXo=zeros(3,3,3,n_modes);
MYo=MXo;
MZo=MXo;
NXo=MXo;
NYo=MXo;
NZo=MXo;
CMXo=zeros(3,3,3,n_modes);
CMYo=CMXo;
CMZo=CMXo;
CNXo=CMXo;
CNYo=CMYo;
CNZo=CMZo;

for jj=1:n_modes
    %extract into cart.
    [Mx(:),My(:),Mz(:)]=xietaphiv2xyzv(isProlate,c,Mfull(:,:,jj),xep(:,1:3));
    MXo(:,:,:,jj)=Mx;
    MYo(:,:,:,jj)=My;
    MZo(:,:,:,jj)=Mz;
    [CMXo(:,:,:,jj),CMYo(:,:,:,jj),CMZo(:,:,:,jj)]=curl(X,Y,Z,Mx,My,Mz);
    [Mx(:),My(:),Mz(:)]=xietaphiv2xyzv(isProlate,c,Nfull(:,:,jj),xep(:,1:3));
    NXo(:,:,:,jj)=Mx;
    NYo(:,:,:,jj)=My;
    NZo(:,:,:,jj)=Mz;
    [CNXo(:,:,:,jj),CNYo(:,:,:,jj),CNZo(:,:,:,jj)]=curl(X,Y,Z,Mx,My,Mz);
end
%% test 2/3
tol=1e-7;
if sqrt(sum(abs(diff([squeeze(CMX(2,2,2,:)),squeeze(NX(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'a: prolate curl M in X components passed.'])
else
    warning(['test ' num2str(testnum) 'a: prolate curl M in X not passed.']);
    warn_state=true;
end


if sqrt(sum(abs(diff([squeeze(CMXo(2,2,2,:)),squeeze(NXo(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'b: oblate curl M in X components passed.'])
else
    warning(['test ' num2str(testnum) 'b: oblate curl M in X not passed']);
    warn_state=true;
end
testnum=testnum+1;


if sqrt(sum(abs(diff([squeeze(CNX(2,2,2,:)),squeeze(MX(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'a: prolate curl N in X components passed.'])
else
    warning(['test va: prolate curl N in X not passed']);
    warn_state=true;
end

if sqrt(sum(abs(diff([squeeze(CNXo(2,2,2,:)),squeeze(MXo(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'b: oblate curl N in X components passed.'])
else
    warning(['test ' num2str(testnum) 'b: oblate curl N in X not passed']);
    warn_state=true;
end
testnum=testnum+1;

%% test 4/5

if sqrt(sum(abs(diff([squeeze(CMY(2,2,2,:)),squeeze(NY(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'a: prolate curl M in Y components passed.'])
else
    warning(['test ' num2str(testnum) 'a: prolate curl M in Y not passed']);
    warn_state=true;
end

if sqrt(sum(abs(diff([squeeze(CMYo(2,2,2,:)),squeeze(NYo(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'b: oblate curl M in Y components passed.'])
else
    warning(['test ' num2str(testnum) 'b: oblate curl M in Y not passed']);
    warn_state=true;
end
testnum=testnum+1;


if sqrt(sum(abs(diff([squeeze(CNY(2,2,2,:)),squeeze(MY(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'a: prolate curl N in Y components passed.'])
else
    warning(['test 5a: prolate curl N in Y not passed']);
    warn_state=true;
end

if sqrt(sum(abs(diff([squeeze(CNYo(2,2,2,:)),squeeze(MYo(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'b: oblate curl N in Y components passed.'])
else
    warning(['test ' num2str(testnum) 'b: oblate curl N in Y not passed']);
    warn_state=true;
end
testnum=testnum+1;

%% test 6a test 6b test 7a test 7b
if sqrt(sum(abs(diff([squeeze(CMZ(2,2,2,:)),squeeze(NZ(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'a: prolate curl M in Z components passed.'])
else
    warning(['test ' num2str(testnum) 'a: prolate curl M in Z not passed']);
    warn_state=true;
end

if sqrt(sum(abs(diff([squeeze(CMZo(2,2,2,:)),squeeze(NZo(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'b: oblate curl M in Z components passed.'])
else
    warning(['test ' num2str(testnum) 'b: oblate curl M in Z not passed']);
    warn_state=true;
end
testnum=testnum+1;


if sqrt(sum(abs(diff([squeeze(CNZ(2,2,2,:)),squeeze(MZ(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'a: prolate curl N in Z components passed.'])
else
    warning(['test ' num2str(testnum) 'a: prolate curl N in Z not passed']);
    warn_state=true;
end

if sqrt(sum(abs(diff([squeeze(CNZo(2,2,2,:)),squeeze(MZo(2,2,2,:))],1,2).^2)))<tol;
    disp(['test ' num2str(testnum) 'b: oblate curl N in Z components passed.'])
else
    warning(['test ' num2str(testnum) 'b: oblate curl N in Z not passed']);
    warn_state=true;
end
testnum=testnum+1;
%%
if warn_state
    error('Spheroidal wavefunction formulation is WRONG!!! Go back and do it again')
end

%% test 8a test 8b test 9a test 9b
disp(['------------BATCH ' num2str(batchno) '-----------'])
batchno=batchno+1;
% first test xietaphi2rtp:

isProlate=true;
c=1.1;
xep=[1.2,.1,.2];

[rtp]=xietaphi2rtp(isProlate,c,xep);
rtpr=[0.7379024325749306,1.390942827002418,0.2000000000000000];

if sqrt(sum(abs(rtp-rtpr).^2))<1e-14;
    disp(['test ' num2str(testnum) 'a: prolate xietaphi2rtp passed.'])
else
    warning(['test ' num2str(testnum) 'a: prolate xietaphi2rtp not passed.']);
    warn_state=true;
end

isProlate=false;
[rtp]=xietaphi2rtp(isProlate,c,xep);
rtpr=[1.714730299493189,1.493740057774045,0.2000000000000000];

if sqrt(sum(abs(rtp-rtpr).^2))<1e-14;
    disp(['test ' num2str(testnum) 'b: oblate xietaphi2rtp passed.'])
else
    warning(['test ' num2str(testnum) 'b: oblate xietaphi2rtp not passed.']);
    warn_state=true;
end
testnum=testnum+1;

% when c is small it should become the VSWF


c=1e-10;
isProlate=true;
rtp=xietaphi2rtp(isProlate,c,[1/c+1,-.1,.4]);
rtpr=[1.000000000100000,1.670963747956456,.4];

if sqrt(sum(abs(rtp-rtpr).^2))<1e-14;
    disp(['test ' num2str(testnum) 'a: prolate xietaphi2rtp passed for c=10^-10.'])
else
    warning(['test ' num2str(testnum) 'a: prolate xietaphi2rtp not passed for c=10^-10.']);
    warn_state=true;
end

isProlate=false;
rtp=xietaphi2rtp(isProlate,c,[1/c,-.1,.4]);
rtpr=[1.000000000000000,1.670963747956456,.4];
if sqrt(sum(abs(rtp-rtpr).^2))<1e-14;
    disp(['test ' num2str(testnum) 'b: oblate xietaphi2rtp passed for c=10^-10.'])
else
    warning(['test ' num2str(testnum) 'b: oblate xietaphi2rtp not passed for c=10^-10.']);
    warn_state=true;
end
testnum=testnum+1;
%% test 10a test 10b

isProlate=true;
c=1e-10;
xi=(abs(randn(10,1))/c+1);
eta=2*(rand(10,1))-1;
phi=((rand(10,1))-.5)*2*pi;

rtp=xietaphi2rtp(isProlate,c,[xi,eta,phi]);
test_for_this=20;
testM=zeros(test_for_this,1);
testN=zeros(test_for_this,1);
for ii=1:test_for_this
    n=randi(12);
    m=randi([0,n]);
    
    [M,N]=spheroidalvwf(isProlate,n,m,c,xi,eta,phi,3);
    [Ms,Ns]=vswf_vector(n,m,rtp(:,1),rtp(:,2),rtp(:,3),3);
    Ms=Ms*sqrt(n*(n+1));
    Ns=Ns*sqrt(n*(n+1));
    
    % Ms(:,3)=(sign(-m))^m*Ms(:,3);
    % Ms(:,2)=(sign(-(m+1)))^(m+1)*Ms(:,2);
    % 
    % Ns(:,[1,3])=(sign(-m))^(m)*Ns(:,[1,3]); %ee fail, oo pass, oe fail, eo pass
    % Ns(:,2)=sign(-(m+1))^(m+1)*Ns(:,2); %ee fail, oo pass, oe fail, eo pass
    M(:,2)=-M(:,2);
    N(:,2)=-N(:,2);
    
    testM(ii)=all(sum(abs(M-Ms).^2,2)<1e-14);
    testN(ii)=all(sum(abs(N-Ns).^2,2)<1e-14);

end

if all(testM==1)
	disp(['test ' num2str(testnum) 'a: prolate M is in limit to M_spherical passed.'])
else
    warning(['test ' num2str(testnum) 'a: prolate M is in limit to M_spherical failed.']);
    warn_state=true;
end

if all(testN==1)
	disp(['test ' num2str(testnum) 'b: prolate N is in limit to N_spherical passed.'])
else
    warning(['test ' num2str(testnum) 'b: prolate N is in limit to N_spherical failed.']);
    warn_state=true;
end
%% test 10c test 10d

isProlate=false;
xi=(abs(randn(10,1))/c);
eta=2*(rand(10,1))-1;
phi=((rand(10,1))-.5)*2*pi;

rtp=xietaphi2rtp(isProlate,c,[xi,eta,phi]);
test_for_this=20;
testM=zeros(test_for_this,1);
testN=zeros(test_for_this,1);
for ii=1:test_for_this
    n=randi(7);
    m=randi([0,n]);
    
    [M,N]=spheroidalvwf(isProlate,n,m,1e-10,xi,eta,phi,3);
    [Ms,Ns]=vswf_vector(n,m,rtp(:,1),rtp(:,2),rtp(:,3),3);
    Ms=Ms*sqrt(n*(n+1));
    Ns=Ns*sqrt(n*(n+1));
    
    % Ms(:,3)=(sign(-m))^m*Ms(:,3);
    % Ms(:,2)=(sign(-(m+1)))^(m+1)*Ms(:,2);
    % 
    % Ns(:,[1,3])=(sign(-m))^(m)*Ns(:,[1,3]); %ee fail, oo pass, oe fail, eo pass
    % Ns(:,2)=sign(-(m+1))^(m+1)*Ns(:,2); %ee fail, oo pass, oe fail, eo pass
    M(:,2)=-M(:,2);
    N(:,2)=-N(:,2);
    
    testM(ii)=all(sum(abs(M-Ms).^2,2)<1e-14);
    testN(ii)=all(sum(abs(N-Ns).^2,2)<1e-14);

end

if all(testM==1)
	disp(['test ' num2str(testnum) 'c: oblate M is in limit to M_spherical passed.'])
else
    warning(['test ' num2str(testnum) 'c: oblate M is in limit to M_spherical failed.']);
    warn_state=true;
end

if all(testN==1)
	disp(['test ' num2str(testnum) 'd: oblate N is in limit to N_spherical passed.'])
else
    warning(['test ' num2str(testnum) 'd: oblate N is in limit to N_spherical failed.']);
    warn_state=true;
end
%%

disp(['------------BATCH ' num2str(batchno) '-----------'])
batchno=batchno+1;
%% test 11a test 11b

%now to match an expansion.
isProlate=true;

[xep]=xyz2xietaphi(isProlate,cpo,xyzp(1,:));

if sqrt(sum(abs(xep(:,1:3)-xepp).^2))<1e-14;
    disp(['test ' num2str(testnum) 'a: prolate xyz2xietaphi passed.']);
else
    warning(['test ' num2str(testnum) 'a: prolate xyz2xietaphi failed.']);
    warn_state=true;
end

isProlate=false;

[xep]=xyz2xietaphi(isProlate,cpo,xyzo(1,:));

if sqrt(sum(abs(xep(:,1:3)-xepo).^2))<1e-14;
    disp(['test ' num2str(testnum) 'b: oblate xyz2xietaphi passed.']);
else
    warning(['test ' num2str(testnum) 'b: oblate xyz2xietaphi failed.']);
    warn_state=true;
end
testnum=testnum+1;

%% compare with VSWF...

htype=3;

VSWFVALUESM=[[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,-0.00790644068018794-0.07880068291302512*1i,-0.03940034145651256+0.00395322034009397*1i];[0,0,0.09699530477132740];[0,0.00790644068018794-0.07880068291302512*1i,0.03940034145651256+0.00395322034009397*1i];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,-0.00401145551294812-0.01978913132189994*1i,-0.009894565660949971+0.002005727756474060*1i];[0,-0.001163821721085791-0.011599397266127666*1i,0.011599397266127666-0.001163821721085791*1i];[0,0,0.02472958229957887];[0,0.001163821721085791-0.011599397266127666*1i,-0.011599397266127666-0.001163821721085791*1i];[0,-0.00401145551294812+0.01978913132189994*1i,-0.009894565660949971-0.002005727756474060*1i];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,-0.000948771824188867-0.003067121378067395*1i,-0.001533560689033697+0.000474385912094434*1i];[0,-0.000601352560880133-0.002966565317652494*1i,0.0007416413294131235-0.0001503381402200333*1i];[0,-0.0000275856770337835-0.0002749366341705890*1i,0.003986581195473540-0.000399992316989861*1i];[0,0,0.0008289512057281660];[0,0.0000275856770337835-0.0002749366341705890*1i,-0.003986581195473540-0.000399992316989861*1i];[0,-0.000601352560880133+0.002966565317652494*1i,0.0007416413294131235+0.0001503381402200333*1i];[0,0.000948771824188867-0.003067121378067395*1i,0.001533560689033697+0.000474385912094434*1i];[0,0,0];[0,0,0];[0,0,0];[0,-0.0001467740676688521-0.0003471533155307066*1i,-0.0001735766577653533+0.0000733870338344261*1i];[0,-0.0001364161341239728-0.0004409962760463006*1i,0];[0,-0.0000283018441767865-0.0001396173805888146*1i,0.0004188521417664439-0.0000849055325303594*1i];[0,9.67682674551017*10^-6+0.00009644549131798479*1i,0.0003086255722175513-0.0000309658455856325*1i];[0,0,-0.0003754072542871961];[0,-9.67682674551017*10^-6+0.00009644549131798479*1i,-0.0003086255722175513-0.0000309658455856325*1i];[0,-0.0000283018441767865+0.0001396173805888146*1i,0.0004188521417664439+0.0000849055325303594*1i];[0,0.0001364161341239728-0.0004409962760463006*1i,0];[0,-0.0001467740676688521+0.0003471533155307066*1i,-0.0001735766577653533-0.0000733870338344261*1i];[0,0,0];[0,-0.00001689503660320748-0.00003092615706005374*1i,-0.00001546307853002687+8.44751830160374*10^-6*1i];[0,-0.00002004396960703978-0.00004740844630115300*1i,-5.926055787644125*10^-6+2.505496200879972*10^-6*1i];[0,-7.76227279597929*10^-6-0.00002509331772715012*1i,0.00003262131304529516-0.00001009095463477308*1i];[0,1.967960053241688*10^-6+9.708251731609476*10^-6*1i,0.00004126006985934027-8.36383022627717*10^-6*1i];[0,1.025050699734409*10^-6+0.000010216315840065928*1i,-0.00001747527709484961+1.75337619691412*10^-6*1i];[0,0,-0.00004870355594125217];[0,-1.025050699734409*10^-6+0.000010216315840065928*1i,0.00001747527709484961+1.75337619691412*10^-6*1i];[0,1.967960053241688*10^-6-9.708251731609476*10^-6*1i,0.00004126006985934027+8.36383022627717*10^-6*1i];[0,7.76227279597929*10^-6-0.00002509331772715012*1i,-0.00003262131304529516-0.00001009095463477308*1i];[0,-0.00002004396960703978+0.00004740844630115300*1i,-5.926055787644125*10^-6-2.505496200879972*10^-6*1i];[0,0.00001689503660320748-0.00003092615706005374*1i,0.00001546307853002687+8.44751830160374*10^-6*1i]];
VSWFVALUESN=[[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0.12407889679680380-0.01244941542283177*1i,0.06265212997883521-0.00628618091688146*1i,-0.01257236183376291-0.12530425995767042*1i];[0.1018186642179194,-0.1542362887534181,0];[-0.12407889679680380-0.01244941542283177*1i,-0.06265212997883521-0.00628618091688146*1i,0.01257236183376291-0.12530425995767042*1i];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0.04673970120979549-0.00947460649190239*1i,0.02539948321329784-0.00514873014406954*1i,-0.01029746028813907-0.05079896642659568*1i];[0.05479294200575087-0.00549763186874412*1i,-0.02977580888751173+0.00298754602080753*1i,-0.00298754602080753-0.02977580888751173*1i];[-0.01297966454328079,-0.06348117057517366,0];[-0.05479294200575087-0.00549763186874412*1i,0.02977580888751173+0.00298754602080753*1i,0.00298754602080753-0.02977580888751173*1i];[0.04673970120979549+0.00947460649190239*1i,0.02539948321329784+0.00514873014406954*1i,-0.01029746028813907+0.05079896642659568*1i];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0,0,0];[0.009658927381442543-0.002987856371427135*1i,0.005386809742329254-0.001666335523052712*1i,-0.003332671046105424-0.010773619484658507*1i];[0.01401338687494316-0.00284065415101649*1i,-0.002605101165649955+0.000528080149792586*1i,-0.002112320599170346-0.010420404662599820*1i];[0.002597477559520569-0.000260617059183813*1i,-0.014003328708102789+0.001405019394032269*1i,-0.0000968978892436048-0.0009657468074553647*1i];[-0.01218242331783837,-0.002911787230113295,0];[-0.002597477559520569-0.000260617059183813*1i,0.014003328708102789+0.001405019394032269*1i,0.0000968978892436048-0.0009657468074553647*1i];[0.01401338687494316+0.00284065415101649*1i,-0.002605101165649955-0.000528080149792586*1i,-0.002112320599170346+0.010420404662599820*1i];[-0.009658927381442543-0.002987856371427135*1i,-0.005386809742329254-0.001666335523052712*1i,0.003332671046105424-0.010773619484658507*1i];[0,0,0];[0,0,0];[0,0,0];[0.001366561773898122-0.000577773050990919*1i,0.0007714780790036655-0.0003261757002078936*1i,-0.000652351400415787-0.001542956158007331*1i];[0.002314630169881401-0.000715999015984398*1i,0,-0.000606314572788575-0.001960050183428981*1i];[0.0010992018036341259-0.0002228192366458705*1i,-0.001861628457861956+0.000377370770797152*1i,-0.0001257902569323841-0.0006205428192873186*1i];[-0.001518622646579025+0.000152370505266045*1i,-0.0013717159130213934+0.0001376306663273959*1i,0.0000430095832273112+0.0004286612228191854*1i];[-0.001458079042687761,0.001668533494711268,0];[0.001518622646579025+0.000152370505266045*1i,0.0013717159130213934+0.0001376306663273959*1i,-0.0000430095832273112+0.0004286612228191854*1i];[0.0010992018036341259+0.0002228192366458705*1i,-0.001861628457861956-0.000377370770797152*1i,-0.0001257902569323841+0.0006205428192873186*1i];[-0.002314630169881401-0.000715999015984398*1i,0,0.000606314572788575-0.001960050183428981*1i];[0.001366561773898122+0.000577773050990919*1i,0.0007714780790036655+0.0003261757002078936*1i,-0.000652351400415787+0.001542956158007331*1i];[0,0,0];[0.0001460882053932764-0.0000798083503431580*1i,0.00008302744145697960-0.00004535809799330751*1i,-0.0000907161959866150-0.0001660548829139592*1i];[0.0002799335376187425-0.0001183540014025884*1i,0.00003181935919318925-0.00001345300929147420*1i,-0.0001076240743317936-0.0002545548735455140*1i];[0.0001975586419722400-0.0000611120493856630*1i,-0.0001751568521015163+0.0000541823637225105*1i,-0.00004167874132500809-0.00013473604007808948*1i];[-0.00011464899035329189+0.00002324050090554925*1i,-0.0002215417860101439+0.0000449087433087710*1i,0.00001056676313147553+0.00005212747906121033*1i];[-0.0002412978831977032+0.0000242105439855549*1i,0.00009383173881705699-9.41457674541705*10^-6*1i,5.50390640501305*10^-6+0.00005485547807766409*1i];[0.00006188873263760771,0.0002615088342083211,0];[0.0002412978831977032+0.0000242105439855549*1i,-0.00009383173881705699-9.41457674541705*10^-6*1i,-5.50390640501305*10^-6+0.00005485547807766409*1i];[-0.00011464899035329189-0.00002324050090554925*1i,-0.0002215417860101439-0.0000449087433087710*1i,0.00001056676313147553-0.00005212747906121033*1i];[-0.0001975586419722400-0.0000611120493856630*1i,0.0001751568521015163+0.0000541823637225105*1i,0.00004167874132500809-0.00013473604007808948*1i];[0.0002799335376187425+0.0001183540014025884*1i,0.00003181935919318925+0.00001345300929147420*1i,-0.0001076240743317936+0.0002545548735455140*1i];[-0.0001460882053932764-0.0000798083503431580*1i,-0.00008302744145697960-0.00004535809799330751*1i,0.0000907161959866150-0.0001660548829139592*1i]];


[nn,mm]=meshgrid([1:5]',[-5:5]');
findx=abs(mm)>nn;
VSWFVALUESM(findx,:)=[];
VSWFVALUESN(findx,:)=[];
nn(findx)=[];
mm(findx)=[];

isProlate=true;
xep=rtp2xietaphi(isProlate,(c),[11/10,pi/3,1/10]);

[M,N]=spheroidalvwf(isProlate,nn,mm,(c),xep(:,1),xep(:,2),xep(:,3),htype);

%note: the theta/phi unit vector is minus between eta and theta coordiante
%systems so...

M(:,2)=-M(:,2);
N(:,2)=-N(:,2);

M_test=(sqrt(nn(:).*(nn(:)+1)).*VSWFVALUESM-M)./abs(VSWFVALUESM);
N_test=(sqrt(nn(:).*(nn(:)+1)).*VSWFVALUESN-N)./abs(VSWFVALUESN);

M_test(isnan(M_test))=0;
N_test(isnan(N_test))=0;
M_test(isinf(M_test))=0;
N_test(isinf(N_test))=0;

if any(abs(M_test(:))>1e-13)
    warning(['test ' num2str(testnum) 'a: M_sph does not appropriately match M for c=1e-10, test failed.']);
else
    disp(['test ' num2str(testnum) 'a: M_sph does appropriately match M for c=1e-10, test passed.']);
    
end
if any(abs(N_test(:))>1e-13)
    warning(['test ' num2str(testnum) 'b: N_sph does not appropriately match N for c=1e-10, test failed.']);
else
    disp(['test ' num2str(testnum) 'b: N_sph does appropriately match N for c=1e-10, test passed.']);
    
end
testnum=testnum+1;

%%
warning('off','MATLAB:rankDeficientMatrix')
htype=3;

% say we want to find the U matrix coefficients. we need to match the orthogonal
% surface for it to yield the U matrix. 

% say we want a minor axis of
minor_axis=100;
% with an aspect ratio of
aspect_ratio=1.5;
isProlate=true;
% when aspect_ratio>1 this corresponds to c of:
c2=sqrt((aspect_ratio.^2-1)*minor_axis);

%xi on this surface is:
xi0=sqrt((minor_axis/c2).^2+1);

xep2=rand(10000,3);
xep2(:,1)=(xep2(:,1)*2-1)/c2/10+xi0;
xep2(:,2)=cos(xep2(:,2)*pi);
xep2(:,3)=pi*(xep2(:,3)*2-1);

xyz2=xietaphi2xyz(isProlate,c2,xep2);

xep=xyz2xietaphi(isProlate,c,xyz2);

m0=randi(7)-1;
n0=m0+randi(2)-1;


u=spheroidal_u_coefficients(isProlate,mod(n0(1)+m0(1),2),m0(1),c2,10);

n=n0+[0:2:2*(size(u,1)-1)];
m=m0*ones(size(n));

[M0,N0]=spheroidalvwf(isProlate,n,m,c,xep(:,1),xep(:,2),xep(:,3),htype);
% we know that this will couple to odd ns of a different c

[M,N]=spheroidalvwf(isProlate,n,m,c2,xep2(:,1),xep2(:,2),xep2(:,3),htype);

UcouplingM0=real(M.'\M0.');
UcouplingN0=real(N.'\N0.');

uT=u(:,1:size(UcouplingM0,2));
% note: sign BS here...
uf=(-1).^([1:size(uT,2)]+[1:size(uT,1)]');

uT=(uf.*uT);

mean_difference_in_absU=abs(mean((uT)-(UcouplingM0)));

if any(mean_difference_in_absU(:,1:7)>15e-5*c2)
    warning(['test ' num2str(testnum) 'a: Prolate |U| matrix for c2=' num2str(c2) ' failed.']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'a: Prolate |U| matrix for c2=' num2str(c2) ' passed.']);
end

M0u=uT*M0;
N0u=uT*N0;


UcouplingMc2=real(M0u.'\M.');
UcouplingNc2=real(N0u.'\N.');

col_error_M=sqrt(sum(abs(UcouplingMc2-eye(size(UcouplingMc2))).^2));
col_error_N=sqrt(sum(abs(UcouplingNc2-eye(size(UcouplingNc2))).^2));

if any(col_error_M(:,1:7)>15e-5*c2)
    warning(['test ' num2str(testnum) 'b: Prolate conversion not equal for, M, c2=' num2str(c2) ', test failed.']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'b: Prolate conversion approximately equal for, M, c2=' num2str(c2) ', test passed.']);
end

if any(col_error_N(:,1:7)>15e-5*c2*10)
    warning(['test ' num2str(testnum) 'c: Prolate conversion not equal for, N, c2=' num2str(c2) ', test failed.']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'c: Prolate conversion approximately equal for, N, c2=' num2str(c2) ', test passed.']);
end

c3=c2/2;

xep3=xyz2xietaphi(isProlate,c3,xyz2);

[M3,N3]=spheroidalvwf(isProlate,n,m,c3,xep3(:,1),xep3(:,2),xep3(:,3),htype);

u=spheroidal_u_coefficients(isProlate,mod(n0(1)+m0(1),2),m0(1),c2,20);
u3=spheroidal_u_coefficients(isProlate,mod(n0(1)+m0(1),2),m0(1),c3,20);

UcouplingM3=real(M.'\M3.');
UcouplingN3=real(N.'\N3.');

ud3=u3*u';

uf=(-1).^([1:size(ud3,2)]+[1:size(ud3,1)]');

ud3=(uf.*ud3);

mean_difference_in_absU=abs(mean((ud3(1:size(UcouplingM3,1),1:size(UcouplingM3,2)))-(UcouplingM3)));

if any(mean_difference_in_absU(:,1:7)>15e-5*c3*c2*10)
    warning(['test ' num2str(testnum) 'd: Prolate conversion not? OK(?) for c2=' num2str(c2) ' to c3=' num2str(c3) ', test failed. NOTE THIS WARNING IS NOT NECESSARILY INDICATIVE OF A PROBLEM AS THIS CAN FAIL FOR REASONS THAT ARE UNRELATED TO IMPLEMENTATION.']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'd: Prolate conversion OK(?) for c2=' num2str(c2) ' to c3=' num2str(c3) ', test passed.']);
end

%% do the same for oblate:

htype=3;

% say we want to find the U matrix coefficients. we need to match the orthogonal
% surface for it to yield the U matrix. 

% say we want a minor axis of
minor_axis=100;
% with an aspect ratio of
aspect_ratio=1/1.5;
isProlate=false;
% when aspect_ratio>1 this corresponds to c of:
c2=sqrt(-(aspect_ratio.^2-1)*minor_axis);

%xi on this surface is:
xi0=sqrt((minor_axis/c2).^2);

xep2=rand(10000,3);
xep2(:,1)=(xep2(:,1)*2-1)/c2/10+xi0;
xep2(:,2)=cos(xep2(:,2)*pi);
xep2(:,3)=pi*(xep2(:,3)*2-1);

xyz2=xietaphi2xyz(isProlate,c2,xep2);

xep=xyz2xietaphi(isProlate,c,xyz2);

m0=randi(7)-1;
n0=m0+randi(2)-1;


u=spheroidal_u_coefficients(isProlate,mod(n0(1)+m0(1),2),m0(1),c2,10);

n=n0+[0:2:2*(size(u,1)-1)];
m=m0*ones(size(n));

[M0,N0]=spheroidalvwf(isProlate,n,m,c,xep(:,1),xep(:,2),xep(:,3),htype);
% we know that this will couple to odd ns of a different c

[M,N]=spheroidalvwf(isProlate,n,m,c2,xep2(:,1),xep2(:,2),xep2(:,3),htype);

UcouplingM0=real(M.'\M0.');
UcouplingN0=real(N.'\N0.');

uT=u(:,1:size(UcouplingM0,2));
% note: sign BS here...
uf=(-1).^([1:size(uT,2)]+[1:size(uT,1)]');

uT=(uf.*uT);

mean_difference_in_absU=abs(mean((uT)-(UcouplingM0)));

if any(mean_difference_in_absU(:,1:7)>15e-5*c2)
    warning(['test ' num2str(testnum) 'a: Oblate |U| matrix for c2=' num2str(c2) ' failed.']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'a: Oblate |U| matrix for c2=' num2str(c2) ' passed.']);
end

M0u=uT*M0;
N0u=uT*N0;


UcouplingMc2=real(M0u.'\M.');
UcouplingNc2=real(N0u.'\N.');

col_error_M=sqrt(sum(abs(UcouplingMc2-eye(size(UcouplingMc2))).^2));
col_error_N=sqrt(sum(abs(UcouplingNc2-eye(size(UcouplingNc2))).^2));

if any(col_error_M(:,1:7)>15e-5*c2)
    warning(['test ' num2str(testnum) 'b: Oblate conversion not equal for, M, c2=' num2str(c2) ', test failed.']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'b: Oblate conversion approximately equal for, M, c2=' num2str(c2) ', test passed.']);
end

if any(col_error_N(:,1:7)>15e-5*c2*10)
    warning(['test ' num2str(testnum) 'c: Oblate conversion not equal for, N, c2=' num2str(c2) ', test failed.']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'c: Oblate conversion approximately equal for, N, c2=' num2str(c2) ', test passed.']);
end

c3=c2/2;

xep3=xyz2xietaphi(isProlate,c3,xyz2);

[M3,N3]=spheroidalvwf(isProlate,n,m,c3,xep3(:,1),xep3(:,2),xep3(:,3),htype);

u=spheroidal_u_coefficients(isProlate,mod(n0(1)+m0(1),2),m0(1),c2,20);
u3=spheroidal_u_coefficients(isProlate,mod(n0(1)+m0(1),2),m0(1),c3,20);

UcouplingM3=real(M.'\M3.');
UcouplingN3=real(N.'\N3.');

ud3=u3*u';

uf=(-1).^([1:size(ud3,2)]+[1:size(ud3,1)]');

ud3=(uf.*ud3);

mean_difference_in_absU=abs(mean((ud3(1:size(UcouplingM3,1),1:size(UcouplingM3,2)))-(UcouplingM3)));

if any(mean_difference_in_absU(:,1:7)>15e-5*c3*c2*10)
    warning(['test ' num2str(testnum) 'd: Oblate conversion not? OK(?) for c2=' num2str(c2) ' to c3=' num2str(c3) ', test failed. NOTE THIS WARNING IS NOT NECESSARILY INDICATIVE OF A PROBLEM AS THIS CAN FAIL FOR REASONS THAT ARE UNRELATED TO IMPLEMENTATION.']);
    % warn_state=true;
else
    disp(['test ' num2str(testnum) 'd: Oblate conversion OK(?) for c2=' num2str(c2) ' to c3=' num2str(c3) ', test passed.']);
end
testnum=testnum+1;

%% test with "x-polarised" plane wave
disp(['------------BATCH ' num2str(batchno) '-----------'])
batchno=batchno+1;



isProlate=true;
cpw=5;

numpts=1000;
xeppw=[1/cpw,0,0].*ones(numpts,1);

xeppw(:,1)=abs(xeppw(:,1)+1*(rand(numpts,1)-.5)/10)+1;
xeppw(:,2)=linspace(-1,1,size(xeppw,1))';
% xeppw(:,3)=rand(numpts,1)*2*pi-pi;

xyz=xietaphi2xyz(isProlate,cpw,xeppw);

k=[0,0,1];
phase_factor=exp(-1i*xyz(:,3));

polarisation=[1,0,0];

Exi=xyzv2xietaphiv(isProlate,cpw, polarisation.*phase_factor,xyz);
Hxi=xyzv2xietaphiv(isProlate,cpw,cross(k.*ones(size(phase_factor)),polarisation.*phase_factor),xyz); %hfield using convention \pm k, \pm \omega.

nmax=20;
[a,b]=bsc_plane(nmax,1,0,0,[1,0]);
ab_vswf=[full(a(find(a))),full(b(find(b)))];
[n,m]=combined_index(reshape(([1:nmax]'.*([2:nmax+1]'))'+[-1,1]',[],1));

[M,N]=spheroidalvwf(isProlate,n,m,cpw,xeppw(:,1),xeppw(:,2),xeppw(:,3),htype);

ab=[M.',N.';1i*N.',1i*M.']\[Exi(:);Hxi(:)];

ab_from_pmp=reshape(ab,[],2);

% convert back to c=0:

[ueven,neven]=spheroidal_u_coefficients(isProlate,false,1,cpw,10);
[uodd,nodd]=spheroidal_u_coefficients(isProlate,true,1,cpw,10);

ueven=ueven';
uodd=uodd';

ueven(:,neven>n(end))=[];
ueven(neven>n(end),:)=[];
neven(neven>n(end))=[];
uodd(:,nodd>n(end))=[];
uodd(nodd>n(end),:)=[];
nodd(nodd>n(end))=[];
ueven=(-1).^([1:size(ueven,2)]+[1:size(ueven,1)]').*ueven;
uodd=(-1).^([1:size(ueven,2)]+[1:size(ueven,1)]').*uodd;

ab_convertedp=ab_from_pmp*0;

ab_convertedp(1:4:end,:)=ueven*ab_from_pmp(1:4:end,:);
ab_convertedp(2:4:end,:)=ueven*ab_from_pmp(2:4:end,:);
ab_convertedp(3:4:end,:)=uodd*ab_from_pmp(3:4:end,:);
ab_convertedp(4:4:end,:)=uodd*ab_from_pmp(4:4:end,:);

to_compare=ab_convertedp.*sqrt(n.*(n+1));

if any(sqrt(sum(abs(to_compare(1:10,:)-ab_vswf(1:10,:)).^2))>1e-6)
    warning(['test ' num2str(testnum) 'a: point-match plane wave, prolate c=' num2str(cpw) ' plane wave fit to vswf check. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'a: point-match plane wave, prolate c=' num2str(cpw) ' plane wave fit to vswf check. passed.']);
end

%%

isProlate=false;
cpw=5;

numpts=1000;
xeppw=[1/cpw,0,0].*ones(numpts,1);

xeppw(:,1)=abs(xeppw(:,1)+1*(rand(numpts,1)-.5)/10)+1;
xeppw(:,2)=linspace(-1,1,size(xeppw,1))';
% xeppw(:,3)=rand(numpts,1)*2*pi-pi;

xyz=xietaphi2xyz(isProlate,cpw,xeppw);

k=[0,0,1];
phase_factor=exp(-1i*xyz(:,3));

polarisation=[1,0,0];

Exi=xyzv2xietaphiv(isProlate,cpw, polarisation.*phase_factor,xyz);
Hxi=xyzv2xietaphiv(isProlate,cpw,cross(k.*ones(size(phase_factor)),polarisation.*phase_factor),xyz); %hfield using convention \pm k, \pm \omega.

nmax=20;
[n,m]=combined_index(reshape(([1:nmax]'.*([2:nmax+1]'))'+[-1,1]',[],1));

[M,N]=spheroidalvwf(isProlate,n,m,cpw,xeppw(:,1),xeppw(:,2),xeppw(:,3),htype);

ab=[M.',N.';1i*N.',1i*M.']\[Exi(:);Hxi(:)];

ab_from_pmo=reshape(ab,[],2);

% convert back to c=0:

[ueven,neven]=spheroidal_u_coefficients(isProlate,false,1,cpw,10);
[uodd,nodd]=spheroidal_u_coefficients(isProlate,true,1,cpw,10);

ueven=ueven';
uodd=uodd';

ueven(:,neven>n(end))=[];
ueven(neven>n(end),:)=[];
neven(neven>n(end))=[];
uodd(:,nodd>n(end))=[];
uodd(nodd>n(end),:)=[];
nodd(nodd>n(end))=[];
ueven=(-1).^([1:size(ueven,2)]+[1:size(ueven,1)]').*ueven;
uodd=(-1).^([1:size(ueven,2)]+[1:size(ueven,1)]').*uodd;

ab_convertedo=ab_from_pmo*0;

ab_convertedo(1:4:end,:)=ueven*ab_from_pmo(1:4:end,:);
ab_convertedo(2:4:end,:)=ueven*ab_from_pmo(2:4:end,:);
ab_convertedo(3:4:end,:)=uodd*ab_from_pmo(3:4:end,:);
ab_convertedo(4:4:end,:)=uodd*ab_from_pmo(4:4:end,:);

to_compare=ab_convertedo.*sqrt(n.*(n+1));

if any(sqrt(sum(abs(to_compare(1:10,:)-ab_vswf(1:10,:)).^2))>1e-6)
    warning(['test ' num2str(testnum) 'b: point-match plane wave, oblate c=' num2str(cpw) ' plane wave fit to vswf check. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'b: point-match plane wave, oblate c=' num2str(cpw) ' plane wave fit to vswf check. passed.']);
end
%% now test far-field plane wave code.
isProlate=true;
[a_s_ff,b_s_ff]=bsc_plane_spheroidal(isProlate,cpw,nmax,0,0,[1,0]);

ad_difference=ab_from_pmp-full([a_s_ff(find(a_s_ff)),b_s_ff(find(b_s_ff))]);

if any(sqrt(sum(abs(ad_difference(1:20,:)).^2))>1e-7)
    warning(['test ' num2str(testnum) 'c: bsc_plane_spheroidal vs point-match plane wave, prolate c=' num2str(cpw) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'c: bsc_plane_spheroidal vs point-match plane wave, prolate c=' num2str(cpw) '. passed. ']);
end


isProlate=false;
[a_s_ff,b_s_ff]=bsc_plane_spheroidal(isProlate,cpw,nmax,0,0,[1,0]);

ad_difference=ab_from_pmo-full([a_s_ff(find(a_s_ff)),b_s_ff(find(b_s_ff))]);

if any(sqrt(sum(abs(ad_difference(1:20,:)).^2))>1e-7)
    warning(['test ' num2str(testnum) 'd: bsc_plane_spheroidal vs point-match plane wave, prolate c=' num2str(cpw) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'd: bsc_plane_spheroidal vs point-match plane wave, oblate c=' num2str(cpw) '. passed. ']);
end

%% test t-matrix method using point matching

nmax=13;
r=.5;
k_medium=2*pi;
k_particle=2*pi*1.2;

ar=1.5;
ac=[ar.^(-1/3),ar.^(2/3)]*r;
isProlate=ar>=1;

polz=[1,0];

[sT,sR,c]=stmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
[a,b]=bsc_plane_spheroidal(isProlate,k_medium*c,nmax,0,0,polz);

ab_ind=find([a;b]);
a_ind=find(a);
b_ind=find(b);

a=a(a_ind);
b=b(b_ind);

[xi0,~,~]=xyz2xietaphi(isProlate,c,[ac(1),0,0]);

[xi,eta,phi]=meshgrid(xi0+[-1e-8,0,1e-8],linspace(-1,1),0);

p_q=sT(ab_ind,ab_ind)*[a;b];
c_d=sR(ab_ind,ab_ind)*[a;b];

a_int=c_d(1:end/2);
b_int=c_d(end/2+1:end);

p_sca=p_q(1:end/2);
q_sca=p_q(end/2+1:end);

[n,m]=combined_index(a_ind-1);

[Mint,Nint]=spheroidalvwf(isProlate,n,m,k_particle*c,xi(:),eta(:),phi(:),3);
[Msca,Nsca]=spheroidalvwf(isProlate,n,m,k_medium*c,xi(:),eta(:),phi(:),1);

xyz=xietaphi2xyz(isProlate,c,[xi(:),eta(:),phi(:)]);

phase_factor=exp(-1i*k_medium*xyz(:,3));
k=[0,0,1];

Einc=xyzv2xietaphiv(isProlate,k_medium*c, [polz,0].*phase_factor,xyz);
Hinc=xyzv2xietaphiv(isProlate,k_medium*c,cross(k.*ones(size(phase_factor)),[polz,0].*phase_factor),xyz); %hfield using convention \pm k, \pm \omega.

Hsca=1i*(p_sca.'*Nsca+q_sca.'*Msca);
Hint=k_particle/k_medium*1i*(a_int.'*Nint+b_int.'*Mint);

Hext=Hinc+reshape(Hsca.',[],3);
Hint=reshape(Hint.',[],3);

% figure(1)
% imagesc(real(reshape(Hint(:,3),size(xi))))

% figure(2)
% imagesc(real(reshape(Hext(:,3),size(xi))))

% figure(3)
% imagesc(real(reshape(Hint(:,3)-Hext(:,3),size(xi))))
% sqrt(sum(abs(Hint(:,3)-Hext(:,3)).^2))
if (sqrt(sum(abs(Hint(:,3)-Hext(:,3)).^2))>1e-2)
    warning(['test ' num2str(testnum) 'e: External H-field did not match Internal H-field, prolate c=' num2str(c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'e: External H-field did match Internal H-field, prolate c=' num2str(c) '. passed. ']);
end




ar=1/1.5;
ac=[ar.^(-1/3),ar.^(2/3)]*r;
isProlate=ar>=1;

polz=[1,0];

[sT,sR,c]=stmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
[a,b]=bsc_plane_spheroidal(isProlate,k_medium*c,nmax,0,0,polz);

ab_ind=find([a;b]);
a_ind=find(a);
b_ind=find(b);

a=a(a_ind);
b=b(b_ind);

[xi0,~,~]=xyz2xietaphi(isProlate,c,[ac(1),0,0]);

[xi,eta,phi]=meshgrid(xi0+[-1e-8,0,1e-8],linspace(-1,1),0);

p_q=sT(ab_ind,ab_ind)*[a;b];
c_d=sR(ab_ind,ab_ind)*[a;b];

a_int=c_d(1:end/2);
b_int=c_d(end/2+1:end);

p_sca=p_q(1:end/2);
q_sca=p_q(end/2+1:end);

[n,m]=combined_index(a_ind-1);

[Mint,Nint]=spheroidalvwf(isProlate,n,m,k_particle*c,xi(:),eta(:),phi(:),3);
[Msca,Nsca]=spheroidalvwf(isProlate,n,m,k_medium*c,xi(:),eta(:),phi(:),1);

xyz=xietaphi2xyz(isProlate,c,[xi(:),eta(:),phi(:)]);

phase_factor=exp(-1i*k_medium*xyz(:,3));
k=[0,0,1];

Einc=xyzv2xietaphiv(isProlate,k_medium*c, [polz,0].*phase_factor,xyz);
Hinc=xyzv2xietaphiv(isProlate,k_medium*c,cross(k.*ones(size(phase_factor)),[polz,0].*phase_factor),xyz); %hfield using convention \pm k, \pm \omega.

Hsca=1i*(p_sca.'*Nsca+q_sca.'*Msca);
Hint=k_particle/k_medium*1i*(a_int.'*Nint+b_int.'*Mint);

Hext=Hinc+reshape(Hsca.',[],3);
Hint=reshape(Hint.',[],3);

% sqrt(sum(abs(Hint(:,3)-Hext(:,3)).^2))
if (sqrt(sum(abs(Hint(:,3)-Hext(:,3)).^2))>1e-2)
    warning(['test ' num2str(testnum) 'f: External H-field did not match Internal H-field, oblate c=' num2str(c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'f: External H-field did match Internal H-field, oblate c=' num2str(c) '. passed. ']);
end
testnum=testnum+1;

warning('on','MATLAB:rankDeficientMatrix')
%% farfield limits

disp(['------------BATCH ' num2str(batchno) '-----------'])
batchno=batchno+1;


k_medium=2*pi;
nmax=10;
polz=[1,0];
c=1.1;
isProlate=true;

[n,m]=combined_index([0:nmax*(nmax+2)]');

[eta,phi,xi]=meshgrid(linspace(-1,1,9),1*linspace(-pi,pi,5),1e7);

[M,N]=spheroidalvwf(isProlate,n,m,k_medium*c,xi(:),eta(:),phi(:),2);
[Mff,Nff]=spheroidalvwf_farfield(isProlate,n,m,k_medium*c,xi(:),eta(:),phi(:),2);


Mcheck=(1-real([M(:)./Mff(:)]));
Mcheck(isnan(Mcheck))=0;
if ~all(Mcheck<1e-10)
    warning(['test ' num2str(testnum) 'a: Far-field match in limit for M, prolate c=' num2str(k_medium*c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'a: Far-field match in limit for M, prolate c=' num2str(k_medium*c) '. passed. ']);
end



Ncheck=(1-real([N(:)./Nff(:)]));
Ncheck(isnan(Ncheck))=0;

if ~all(Ncheck<1e-10)
    warning(['test ' num2str(testnum) 'b: Far-field match in limit for N, prolate c=' num2str(k_medium*c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'b: Far-field match in limit for N, prolate c=' num2str(k_medium*c) '. passed. ']);
end


isProlate=false;

[eta,phi,xi]=meshgrid(linspace(-1,1,9),linspace(-pi,pi,5),1e7);

[M,N]=spheroidalvwf(isProlate,n,m,k_medium*c,xi(:),eta(:),phi(:),2);
[Mff,Nff]=spheroidalvwf_farfield(isProlate,n,m,k_medium*c,xi(:),eta(:),phi(:),2);

Mcheck=(1-real([M(:)./Mff(:)]));
Mcheck(isnan(Mcheck))=0;
if ~all(Mcheck<1e-10)
    warning(['test ' num2str(testnum) 'c: Far-field match in limit for M, oblate c=' num2str(k_medium*c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'c: Far-field match in limit for M, oblate c=' num2str(k_medium*c) '. passed. ']);
end

Ncheck=(1-real([N(:)./Nff(:)]));
Ncheck(isnan(Ncheck))=0;
if ~all(Ncheck<1e-10)
    warning(['test ' num2str(testnum) 'd: Far-field match in limit for N, oblate c=' num2str(k_medium*c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'd: Far-field match in limit for N, oblate c=' num2str(k_medium*c) '. passed. ']);
end
testnum=testnum+1;


%% EBCM check


disp(['------------BATCH ' num2str(batchno) '-----------'])
batchno=batchno+1;

warning('off','MATLAB:rankDeficientMatrix');

nmax=13;
r=.5;
k_medium=2*pi;
k_particle=2*pi*1.2;

ar=1.5;
ac=[ar.^(-1/3),ar.^(2/3)]*r;

[sT,sR,c]=stmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
[sTe,sRe,ce]=stmatrix_spheroid_ebcm(nmax,k_medium,k_particle,ac);

[n,m]=combined_index([0:nmax*(nmax+2)]');

mask=ones(size(sT));

meq0=find(m==0);

mask(meq0,meq0)=0;
mask(nmax*(nmax+2)+1+meq0,meq0)=0;
mask(meq0,nmax*(nmax+2)+1+meq0)=0;
mask(nmax*(nmax+2)+1+meq0,nmax*(nmax+2)+1+meq0)=0;


Tprolcheck=sqrt(sum(abs(mask.*(sT-sTe)).^2,'all'));
if Tprolcheck>1e-5
    warning(['test ' num2str(testnum) 'a: EBCM match for prolate, c=' num2str(k_medium*c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'a: EBCM match for prolate, c=' num2str(k_medium*c) '. passed. ']);
end


ar=1/1.5;
ac=[ar.^(-1/3),ar.^(2/3)]*r;


[sT,sR,c]=stmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
[sTe,sRe,ce]=stmatrix_spheroid_ebcm(nmax,k_medium,k_particle,ac);

[n,m]=combined_index([0:nmax*(nmax+2)]');

mask=ones(size(sT));

meq0=find(m==0);

mask(meq0,meq0)=0;
mask(nmax*(nmax+2)+1+meq0,meq0)=0;
mask(meq0,nmax*(nmax+2)+1+meq0)=0;
mask(nmax*(nmax+2)+1+meq0,nmax*(nmax+2)+1+meq0)=0;


Toblcheck=sqrt(sum(abs(mask.*(sT-sTe)).^2,'all'));
if Toblcheck>3e-5
    warning(['test ' num2str(testnum) 'b: EBCM match for oblate, c=' num2str(k_medium*c) '. failed. ']);
    warn_state=true;
else
    disp(['test ' num2str(testnum) 'b: EBCM match for oblate, c=' num2str(k_medium*c) '. passed. ']);
end
warning('on','MATLAB:rankDeficientMatrix');


%% end test

if warn_state
    warning('TEST FINISHED WITH NOTABLE (UNEXPECTED) WARNINGS.')
else
    disp('TEST FINISHED WITHOUT NOTABLE (UNEXPECTED) WARNINGS.')
end