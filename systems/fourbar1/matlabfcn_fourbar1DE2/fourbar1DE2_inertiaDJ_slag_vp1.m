% Calculate time derivative of joint inertia matrix for
% fourbar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% MqD [1x1]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:05
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1DE2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1DE2_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1DE2_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:05:04
% EndTime: 2020-04-24 20:05:07
% DurationCPUTime: 0.76s
% Computational Cost: add. (1889->73), mult. (2709->160), div. (84->11), fcn. (691->4), ass. (0->75)
t69 = pkin(2) * qJD(1);
t93 = pkin(4) ^ 2;
t48 = pkin(2) ^ 2;
t49 = pkin(1) ^ 2;
t42 = cos(qJ(1));
t81 = pkin(2) * t42;
t71 = -0.2e1 * pkin(1) * t81 + t49;
t36 = t48 + t71;
t32 = 0.1e1 / t36;
t92 = 0.4e1 * t32;
t41 = sin(qJ(1));
t66 = t41 * t69;
t62 = pkin(1) * t66;
t70 = pkin(3) ^ 2 - t93;
t31 = t36 - t70;
t85 = -pkin(3) - pkin(4);
t20 = (pkin(2) - t85) * (pkin(2) + t85) + t71;
t84 = -pkin(3) + pkin(4);
t21 = (pkin(2) - t84) * (pkin(2) + t84) + t71;
t73 = t20 * t21;
t50 = sqrt(-t73);
t63 = -pkin(1) + t81;
t80 = t41 * pkin(2);
t11 = t31 * t80 - t50 * t63;
t15 = 0.1e1 / t50;
t82 = pkin(1) * t32;
t67 = t15 * t82;
t3 = -t11 * t67 + 0.1e1;
t46 = 0.1e1 / pkin(3);
t91 = t3 * t46;
t18 = 0.1e1 / t21;
t90 = 0.1e1 / t20 ^ 2 * t18;
t52 = t36 ^ 2;
t33 = 0.1e1 / t52;
t89 = t33 * Icges(4,3);
t30 = t36 + t70;
t38 = pkin(1) * t42 - pkin(2);
t12 = -pkin(1) * t30 * t41 + t38 * t50;
t10 = t12 ^ 2;
t61 = t48 * t62;
t88 = t10 * t61 * t89;
t34 = 0.1e1 / t52 ^ 2;
t87 = -t32 / 0.2e1;
t86 = t32 / 0.2e1;
t76 = rSges(4,2) * t41;
t26 = -rSges(4,1) * t63 + pkin(2) * t76;
t78 = rSges(4,1) * t41;
t27 = rSges(4,2) * t63 + pkin(2) * t78;
t6 = t26 * t50 + t27 * t31;
t7 = -t26 * t31 + t27 * t50;
t83 = t6 ^ 2 + t7 ^ 2;
t79 = rSges(3,1) * t41;
t77 = rSges(3,2) * t41;
t13 = (-t20 - t21) * t62;
t75 = t13 * t15;
t16 = 0.1e1 / t20;
t74 = t16 * t18;
t72 = t41 * t50;
t28 = rSges(3,2) * t63 + pkin(2) * t79;
t29 = -rSges(3,1) * t63 + pkin(2) * t77;
t8 = t28 * t50 + t29 * t30;
t65 = t8 * t86;
t9 = -t28 * t30 + t29 * t50;
t64 = t9 * t87;
t40 = t41 ^ 2;
t60 = t12 * (t38 * t75 + (-0.2e1 * t49 * t40 * pkin(2) + (-t30 * t42 - t72) * pkin(1)) * qJD(1)) * t48 * t74;
t59 = t83 * t61;
t58 = t34 * t59;
t25 = (rSges(4,2) * t42 + t78) * t69;
t24 = (rSges(3,2) * t42 + t79) * t69;
t23 = (rSges(3,1) * t42 - t77) * t69;
t22 = (rSges(4,1) * t42 - t76) * t69;
t19 = 0.1e1 / t21 ^ 2;
t1 = -(-t63 * t75 + (0.2e1 * t48 * t40 * pkin(1) + (t31 * t42 + t72) * pkin(2)) * qJD(1)) * t67 + (-t13 / t73 * t82 + 0.2e1 * t49 * t33 * t66) * t11 * t15;
t2 = [-0.2e1 * t60 * t89 + t74 * t88 * t92 + 0.2e1 * t3 * Icges(3,3) * t1 + 0.2e1 * m(3) * ((t64 * t91 - t80) * (-t42 * t69 + (t1 * t64 + ((-t30 * t23 + t24 * t50 + t29 * t75) * t87 + (t28 * t32 + t33 * t9) * t62) * t3) * t46) + (t65 * t91 + t81) * (-t66 + (t1 * t65 + ((t23 * t50 + t30 * t24 + t28 * t75) * t86 + (t29 * t32 - t33 * t8) * t62) * t3) * t46)) + m(4) * (-t83 * t34 * t60 + (t58 * t90 + (t19 * t58 + (t59 * t92 + (-t6 * (t31 * t22 + t25 * t50 + t26 * t75 + 0.2e1 * t27 * t62) - t7 * (t22 * t50 - t31 * t25 - 0.2e1 * t26 * t62 + t27 * t75)) * t48) * t18 * t34) * t16) * t10) / t93 / 0.2e1 + 0.2e1 * (t19 * t16 + t90) * t88;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t2(1);];
Mq = res;
