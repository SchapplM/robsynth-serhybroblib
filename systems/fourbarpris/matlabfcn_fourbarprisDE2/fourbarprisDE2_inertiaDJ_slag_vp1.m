% Calculate time derivative of joint inertia matrix for
% fourbarprisDE2
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
% MqD [1x1]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisDE2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:45:03
% EndTime: 2020-05-07 09:45:07
% DurationCPUTime: 0.48s
% Computational Cost: add. (1030->81), mult. (1356->123), div. (30->11), fcn. (42->2), ass. (0->59)
t28 = (qJ(1) + pkin(3));
t41 = pkin(2) ^ 2;
t43 = (pkin(1) ^ 2);
t93 = (t41 ^ 2 + t43 ^ 2);
t92 = -2 * t43;
t25 = -2 * qJ(1) * rSges(3,3);
t32 = rSges(3,3) ^ 2;
t36 = (rSges(3,1) ^ 2);
t70 = -t32 - t36;
t61 = t25 + t70;
t64 = -pkin(2) - t28;
t19 = pkin(1) - t64;
t14 = 0.1e1 / t19;
t63 = -pkin(2) + t28;
t20 = pkin(1) - t63;
t15 = 0.1e1 / t20;
t91 = t14 * t15;
t21 = pkin(1) + t63;
t16 = 0.1e1 / t21;
t22 = pkin(1) + t64;
t17 = 0.1e1 / t22;
t90 = t16 * t17;
t72 = t21 * t22;
t62 = t20 * t72;
t44 = sqrt(-t19 * t62);
t75 = rSges(3,1) / t44 * (-t62 + (t72 + (t21 - t22) * t20) * t19) * qJD(1) / 0.2e1;
t89 = (-2 * pkin(3) - 4 * qJ(1)) * rSges(3,3) * qJD(1) + t75;
t39 = pkin(3) ^ 2;
t88 = -3 * t39 + t70;
t56 = -t32 / 0.2e1 - t36 / 0.2e1;
t29 = m(2) * rSges(2,2) ^ 2;
t30 = m(2) * rSges(2,1) ^ 2;
t55 = (-t29 - t30 - Icges(3,2) - Icges(2,3));
t74 = pkin(3) * qJ(1);
t86 = (t39 + 2 * t74 + t61) * m(3) + t55;
t85 = -2 * Icges(4,3) - 2 * (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4);
t82 = 2 * t28;
t49 = t28 ^ 2;
t50 = t28 * t49;
t80 = 2 * t50;
t79 = -0.2e1 * qJD(1);
t78 = t55 / 0.2e1 + t85;
t76 = 0.3e1 / 0.2e1 * t39;
t73 = rSges(3,1) * t44;
t31 = qJ(1) ^ 2;
t71 = t31 * rSges(3,3);
t68 = qJD(1) * m(3);
t66 = t28 * t92;
t59 = -(2 * t71) + t73;
t58 = 3 * t39 + t70;
t51 = t55 * t28;
t38 = pkin(3) * t39;
t27 = 1 / t49;
t11 = (pkin(3) - rSges(3,3)) * t68;
t5 = (t38 + ((2 * t31 + t61) * pkin(3)) + (t58 * qJ(1)) + t59) * m(3) + t51;
t4 = (-t38 + ((-4 * t31 + t61) * pkin(3)) + t59 + ((-2 * t31 + t88) * qJ(1))) * m(3) + t51;
t2 = (-t29 / 0.2e1 - t30 / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(2,3) / 0.2e1 + t85) * pkin(3) + qJ(1) * t78 + (t38 / 0.2e1 + (t31 + t56) * pkin(3) - t71 + t73 + (-(pkin(3) * rSges(3,3)) + t56 + t76) * qJ(1)) * m(3);
t1 = (t4 * t82 + (t86 * t92)) * t41 + t5 * t66 + t2 * t80 + (t93 * t86);
t3 = [(((-0.4e1 * t11 * t43 + t89 * m(3) * t82 + 0.2e1 * (t4 + ((t55 + (-8 * t74 - 6 * t31 + t88) * m(3)) * t28)) * qJD(1)) * t41 + t5 * t43 * t79 + (qJD(1) * t55 + ((t58 + 4 * t74) * qJD(1) + t89) * m(3)) * t66 + 0.6e1 * t49 * t2 * qJD(1) + ((t75 + (t56 + t25) * qJD(1)) * m(3) + qJD(1) * t78 + (t76 + ((-rSges(3,3) + 2 * qJ(1)) * pkin(3))) * t68) * t80 + 0.2e1 * t93 * t11) * t27 + t1 / t50 * t79) * t90 * t91 + ((0.1e1 / t22 ^ 2 * t16 - t17 / t21 ^ 2) * t91 + (0.1e1 / t20 ^ 2 * t14 - t15 / t19 ^ 2) * t90) * qJD(1) * t1 * t27;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t3(1);];
Mq = res;
