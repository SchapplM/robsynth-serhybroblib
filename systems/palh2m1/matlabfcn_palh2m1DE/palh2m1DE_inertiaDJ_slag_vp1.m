% Calculate time derivative of joint inertia matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1DE_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1DE_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:13
% DurationCPUTime: 1.51s
% Computational Cost: add. (296->90), mult. (547->157), div. (0->0), fcn. (302->18), ass. (0->55)
t95 = Icges(3,1) - Icges(3,2);
t94 = Icges(4,1) - Icges(4,2);
t36 = sin(qJ(1));
t38 = cos(qJ(3));
t18 = pkin(3) * t38 + pkin(2);
t34 = sin(qJ(3));
t35 = sin(qJ(2));
t39 = cos(qJ(2));
t48 = t34 * t39 + t35 * t38;
t57 = qJD(2) * t39;
t58 = qJD(2) * t35;
t87 = pkin(3) * (qJD(3) * t48 + t34 * t57) + t18 * t58;
t92 = t87 * t36;
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t13 = rSges(6,1) * t33 + rSges(6,2) * t37;
t29 = qJD(2) + qJD(3);
t88 = t13 * t29;
t69 = t34 * t35;
t46 = -pkin(3) * t69 + t18 * t39 + pkin(1);
t5 = rSges(5,1) + t46;
t40 = cos(qJ(1));
t14 = rSges(6,1) * t40 - t36 * rSges(6,2);
t15 = rSges(6,1) * t36 + rSges(6,2) * t40;
t86 = -t14 * t37 + t15 * t33;
t80 = m(6) * pkin(3);
t27 = qJD(4) + qJD(2);
t79 = t27 / 0.2e1;
t78 = rSges(3,3) * m(3);
t77 = rSges(4,3) * m(4);
t76 = rSges(5,3) * m(5);
t16 = rSges(6,1) * t37 - rSges(6,2) * t33;
t7 = pkin(4) + t46;
t75 = t16 + t7;
t32 = qJ(2) + qJ(3);
t23 = sin(t32);
t71 = t23 * t29;
t24 = cos(t32);
t70 = t24 * t29;
t65 = t38 * t39;
t64 = t40 * t13;
t59 = qJD(4) * t7;
t31 = -qJ(4) + qJ(2);
t30 = qJ(4) + qJ(2);
t28 = -qJD(4) + qJD(2);
t56 = t80 / 0.2e1;
t52 = (qJD(3) + t27) * t56;
t22 = qJD(3) + t28;
t25 = qJ(3) + t30;
t26 = qJ(3) + t31;
t43 = -(rSges(4,2) * t77 - Icges(4,6)) * t71 + (rSges(4,1) * t77 + pkin(3) * t76 - Icges(4,5)) * t70 + (t22 * cos(t26) * t56 + cos(t25) * t52) * rSges(6,2) + (sin(t25) * t52 - t22 * sin(t26) * t80 / 0.2e1) * rSges(6,1);
t12 = t16 * qJD(4);
t3 = (-(m(4) * rSges(4,1) + pkin(3) * (m(5) + m(6))) * t34 - m(4) * t38 * rSges(4,2)) * qJD(3) * pkin(2);
t1 = t87 * t40;
t2 = [0.2e1 * m(6) * ((-t14 * t33 - t15 * t37 - t36 * t7) * (t92 + t86 * qJD(4) + (-t7 * t40 + t86) * qJD(1)) + (-t36 * t13 + t40 * t75) * (-t36 * t12 - t1 - qJD(4) * t64 + (-t36 * t75 - t64) * qJD(1))) + 0.2e1 * m(5) * ((-rSges(5,3) * t40 - t36 * t5) * t92 - (-t36 * rSges(5,3) + t40 * t5) * t1) + (-0.2e1 * Icges(4,4) * t23 + t94 * t24) * t71 + (0.2e1 * Icges(4,4) * t24 + t94 * t23) * t70 + (-0.2e1 * Icges(3,4) * t35 + t95 * t39) * t58 + (0.2e1 * Icges(3,4) * t39 + t95 * t35) * t57 + (m(3) * (rSges(3,1) * t39 - rSges(3,2) * t35 + pkin(1)) * (rSges(3,1) * t35 + rSges(3,2) * t39) * qJD(2) + m(4) * (rSges(4,1) * t24 - rSges(4,2) * t23 + pkin(2) * t39 + pkin(1)) * ((rSges(4,1) * t23 + rSges(4,2) * t24) * t29 + pkin(2) * t58)) * (-0.2e1 * t36 ^ 2 - 0.2e1 * t40 ^ 2); (((t76 + t77) * pkin(2) + rSges(3,1) * t78 - Icges(3,5)) * t39 - t35 * (rSges(3,2) * t78 - Icges(3,6))) * qJD(2) + ((t28 * cos(t31) / 0.2e1 + cos(t30) * t79) * rSges(6,2) + (-t28 * sin(t31) / 0.2e1 + sin(t30) * t79) * rSges(6,1)) * pkin(2) * m(6) + t43; 0.2e1 * t3; t43; t3; 0; ((t33 * t87 - t37 * t59) * rSges(6,2) + (-t33 * t59 - t37 * t87) * rSges(6,1)) * m(6); ((t35 * t12 + t13 * t57) * t18 + (qJD(3) * t13 * t65 + (t39 * t12 - t35 * t88) * t34) * pkin(3)) * m(6); (t48 * t12 + (t65 - t69) * t88) * t80; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
