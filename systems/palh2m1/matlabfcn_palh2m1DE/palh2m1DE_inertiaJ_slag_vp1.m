% Calculate joint inertia matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1DE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1DE_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:05
% EndTime: 2020-05-02 23:52:05
% DurationCPUTime: 0.19s
% Computational Cost: add. (167->86), mult. (266->124), div. (0->0), fcn. (120->18), ass. (0->47)
t38 = cos(qJ(3));
t20 = t38 * pkin(3) + pkin(2);
t34 = sin(qJ(3));
t35 = sin(qJ(2));
t39 = cos(qJ(2));
t45 = -t35 * t34 * pkin(3) + t20 * t39 + pkin(1);
t60 = rSges(5,1) + t45;
t59 = m(6) * pkin(3);
t41 = m(5) + m(6);
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t16 = t33 * rSges(6,1) + t37 * rSges(6,2);
t58 = m(6) * t16;
t57 = rSges(3,3) * m(3);
t56 = rSges(4,3) * m(4);
t55 = rSges(5,3) * m(5);
t54 = (m(4) * rSges(4,1) + pkin(3) * t41) * t38;
t53 = t39 * t34;
t52 = (rSges(6,1) ^ 2 + rSges(6,2) ^ 2) * m(6) + Icges(6,3);
t31 = -qJ(4) + qJ(2);
t30 = qJ(4) + qJ(2);
t51 = m(4) * t34 * rSges(4,2);
t50 = t59 / 0.2e1;
t49 = t41 * pkin(3) ^ 2 + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3);
t46 = rSges(6,1) * t37 - rSges(6,2) * t33;
t44 = rSges(3,1) * t39 - rSges(3,2) * t35 + pkin(1);
t32 = qJ(2) + qJ(3);
t25 = sin(t32);
t26 = cos(t32);
t43 = rSges(4,1) * t26 - rSges(4,2) * t25 + t39 * pkin(2) + pkin(1);
t27 = qJ(3) + t30;
t28 = qJ(3) + t31;
t42 = t26 * (rSges(4,2) * t56 - Icges(4,6)) + (rSges(4,1) * t56 + pkin(3) * t55 - Icges(4,5)) * t25 + (sin(t27) + sin(t28)) * rSges(6,2) * t50 + (-cos(t27) * t59 / 0.2e1 + cos(t28) * t50) * rSges(6,1);
t40 = cos(qJ(1));
t36 = sin(qJ(1));
t17 = t40 * rSges(2,1) - t36 * rSges(2,2);
t15 = -t36 * rSges(2,1) - t40 * rSges(2,2);
t9 = pkin(4) + t45;
t8 = -t36 * rSges(3,3) + t44 * t40;
t7 = -t40 * rSges(3,3) - t44 * t36;
t6 = -t40 * rSges(5,3) - t60 * t36;
t5 = -t36 * rSges(5,3) + t60 * t40;
t4 = -t36 * rSges(4,3) + t43 * t40;
t3 = -t40 * rSges(4,3) - t43 * t36;
t2 = -t36 * t16 + (t46 + t9) * t40;
t1 = -t9 * t36 - (rSges(6,1) * t36 + t40 * rSges(6,2)) * t37 - t33 * (t40 * rSges(6,1) - t36 * rSges(6,2));
t10 = [t39 ^ 2 * Icges(3,2) + t26 ^ 2 * Icges(4,2) + Icges(5,2) + Icges(2,3) + Icges(6,3) + (Icges(3,1) * t35 + 0.2e1 * Icges(3,4) * t39) * t35 + (Icges(4,1) * t25 + 0.2e1 * Icges(4,4) * t26) * t25 + m(6) * (t1 ^ 2 + t2 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(3) * (t7 ^ 2 + t8 ^ 2) + m(2) * (t15 ^ 2 + t17 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2); ((t55 + t56) * pkin(2) + rSges(3,1) * t57 - Icges(3,5)) * t35 + t39 * (rSges(3,2) * t57 - Icges(3,6)) + ((sin(t31) / 0.2e1 + sin(t30) / 0.2e1) * rSges(6,2) + (cos(t31) / 0.2e1 - cos(t30) / 0.2e1) * rSges(6,1)) * pkin(2) * m(6) + t42; (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3) + Icges(3,3) + (-0.2e1 * t51 + 0.2e1 * t54 + (m(4) + t41) * pkin(2)) * pkin(2) + t49; t42; (-t51 + t54) * pkin(2) + t49; t49; t46 * t9 * m(6) + t52; (pkin(3) * t53 + t35 * t20) * t58; pkin(3) * (t35 * t38 + t53) * t58; t52;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
