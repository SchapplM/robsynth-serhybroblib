% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [1x1]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbarprisTE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_coriolismatJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_coriolismatJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_coriolismatJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_coriolismatJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisTE_coriolismatJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisTE_coriolismatJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:21
% EndTime: 2020-05-07 09:01:24
% DurationCPUTime: 0.42s
% Computational Cost: add. (1030->91), mult. (1319->123), div. (30->11), fcn. (42->2), ass. (0->63)
t76 = (qJ(1) * rSges(3,3));
t24 = -2 * t76;
t40 = rSges(3,3) ^ 2;
t44 = rSges(3,1) ^ 2;
t78 = (-t40 - t44);
t68 = t24 + t78;
t47 = (pkin(3) ^ 2);
t92 = t78 - 3 * t47;
t91 = ((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4));
t49 = (pkin(2) ^ 2);
t51 = (pkin(1) ^ 2);
t90 = (t49 ^ 2 + t51 ^ 2);
t31 = -t40 / 0.2e1;
t89 = t31 - t44 / 0.2e1;
t32 = (qJ(1) + pkin(3));
t87 = 2 * t32;
t57 = t32 ^ 2;
t58 = t32 * t57;
t85 = 2 * t58;
t84 = 2 * qJ(1);
t83 = 0.3e1 / 0.2e1 * t47;
t38 = qJ(1) ^ 2;
t82 = rSges(3,3) * t38;
t71 = -pkin(2) - t32;
t21 = pkin(1) + t71;
t17 = 1 / t21;
t28 = 1 / t57;
t81 = t17 * t28;
t70 = -pkin(2) + t32;
t20 = pkin(1) + t70;
t80 = t20 * t21;
t18 = pkin(1) - t71;
t19 = pkin(1) - t70;
t69 = t19 * t80;
t52 = sqrt(-(t18 * t69));
t79 = t52 * rSges(3,1);
t77 = -t40 + (3 * t47);
t35 = rSges(2,2) ^ 2 * m(2);
t36 = m(2) * rSges(2,1) ^ 2;
t64 = (-t35 - t36 - Icges(3,2) - Icges(2,3));
t75 = qJ(1) * t64;
t74 = -2 * t32 * t51;
t73 = 2 * t49;
t72 = 2 * t51;
t3 = 0.1e1 / t52 * (-t69 + (t80 + (t20 - t21) * t19) * t18) / 0.2e1;
t67 = -t81 / 0.2e1;
t66 = -(2 * t82) + t79;
t16 = 1 / t20;
t65 = (t16 * t81) / 0.2e1;
t63 = pkin(3) * t84 + t47 + t68;
t46 = pkin(3) * t47;
t62 = -t46 + t66 + ((-2 * t38 + t92) * qJ(1));
t61 = ((-4 * t38 + t68) * m(3) + t64) * pkin(3) + t75;
t60 = t32 * t64;
t59 = -t35 / 0.2e1 - t36 / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(2,3) / 0.2e1 - (2 * Icges(4,3)) - (2 * t91);
t39 = -2 * rSges(3,3);
t25 = -4 * t76;
t15 = 1 / t19;
t14 = 1 / t18;
t13 = (pkin(3) - rSges(3,3)) * m(3);
t4 = (t64 * pkin(3)) + t75 + (t46 + ((2 * t38 + t68) * pkin(3)) + (t77 - t44) * qJ(1) + t66) * m(3);
t2 = t59 * pkin(3) - (qJ(1) * (4 * Icges(4,3) + 4 * t91 - t64)) / 0.2e1 + (t46 / 0.2e1 + (t38 + t89) * pkin(3) - t82 + t79 + (-(rSges(3,3) * pkin(3)) + t83 + t89) * qJ(1)) * m(3);
t1 = [((((t60 + t61) * t73) + 0.6e1 * t57 * t2 + t59 * t85 + (-(2 * t13 * t49) - t4 - t60) * t72 + (2 * t90 * t13) + ((t32 * (((t39 - 8 * qJ(1)) * pkin(3)) + t25 - (6 * t38) + rSges(3,1) * t3 + t92) + t62) * t73 + (((t39 + 4 * qJ(1)) * pkin(3)) + t25 + (-rSges(3,1) + t3) * rSges(3,1) + t77) * t74 + (t83 + ((-rSges(3,3) + t84) * pkin(3)) + t31 + t24 + (-rSges(3,1) / 0.2e1 + t3) * rSges(3,1)) * t85) * m(3)) * t15 * t14 * t65 + (t16 * t15 / (t18 ^ 2) * t67 + (0.1e1 / (t19 ^ 2) * t65 + (0.1e1 / (t20 ^ 2) * t67 + (-(1 / t58 * t17) + (t28 / t21 ^ 2) / 0.2e1) * t16) * t15) * t14) * ((-(t64 * t72) + (t61 * t87) + (t62 * t87 - (t63 * t72)) * m(3)) * t49 + t4 * t74 + t2 * t85 + (t90 * (t63 * m(3) + t64)))) * qJD(1);];
Cq = t1;
