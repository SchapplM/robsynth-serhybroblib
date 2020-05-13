% Calculate joint inertia matrix for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2DE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2DE_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:27
% EndTime: 2020-05-03 01:06:27
% DurationCPUTime: 0.38s
% Computational Cost: add. (135->94), mult. (223->113), div. (0->0), fcn. (42->6), ass. (0->45)
t14 = (m(5) * rSges(5,1));
t58 = 2 * pkin(5) * m(6) + 2 * t14;
t55 = 2 * pkin(1);
t15 = m(6) + m(7);
t57 = t15 * pkin(5) ^ 2;
t16 = m(5) + m(6);
t7 = m(4) + t16;
t56 = (m(7) + t7) * pkin(4) ^ 2;
t11 = cos(qJ(4));
t8 = sin(qJ(4));
t30 = rSges(7,1) * t11 - rSges(7,2) * t8;
t9 = sin(qJ(3));
t54 = pkin(5) * t9;
t53 = (t8 * rSges(7,1) + t11 * rSges(7,2)) * m(7);
t10 = sin(qJ(2));
t50 = pkin(4) * t10;
t49 = m(4) * rSges(4,1);
t48 = m(6) * rSges(6,1);
t47 = rSges(5,2) * m(5);
t45 = rSges(3,3) * m(3);
t44 = rSges(5,3) * m(5);
t12 = cos(qJ(3));
t43 = t12 * pkin(5);
t40 = rSges(7,1) ^ 2 + rSges(7,2) ^ 2;
t39 = t40 * m(7) + Icges(7,3);
t38 = m(3) * pkin(1);
t37 = t16 * pkin(2);
t36 = 2 * t48;
t35 = 2 * m(7);
t34 = t9 * t47;
t33 = pkin(1) + pkin(2);
t32 = pkin(2) + pkin(3);
t31 = -0.2e1 * t34;
t6 = pkin(1) + t32;
t29 = t30 + t6;
t28 = -(rSges(6,3) * m(6)) + t53;
t25 = (rSges(3,1) ^ 2);
t24 = (rSges(5,1) ^ 2);
t22 = (rSges(3,2) ^ 2);
t21 = rSges(5,2) ^ 2;
t19 = pkin(1) ^ 2;
t18 = pkin(2) ^ 2;
t13 = cos(qJ(2));
t3 = t15 * pkin(5) + t14;
t1 = [t33 * t31 - 0.2e1 * t10 * rSges(3,2) * t38 + ((m(3) + t7) * t19) + ((t37 + t48 + t49) * t55) + (t16 * t18) + (pkin(2) * t36) + ((rSges(5,3) ^ 2 + t21) * m(5)) + ((rSges(6,1) ^ 2 + rSges(6,3) ^ 2) * m(6)) + ((rSges(3,3) ^ 2 + t22) * m(3)) + ((rSges(4,1) ^ 2 + rSges(4,3) ^ 2) * m(4)) + Icges(3,1) + Icges(5,1) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(4,2) + Icges(6,2) + Icges(2,3) + Icges(7,3) + (0.2e1 * (-rSges(5,1) * t47 + Icges(5,4)) * t9 + (t33 * t58) + (t29 * t35 + t36) * pkin(5) + ((-t21 + t24) * m(5) - Icges(5,1) + Icges(5,2) + t57) * t12) * t12 + (0.2e1 * (-m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4)) * t10 + 0.2e1 * rSges(3,1) * t38 + (t12 * t58 + t31 + (t7 * t55) + (2 * t37) + (2 * t49) + t36 + (t29 + t43) * t35) * pkin(4) + ((-t22 + t25) * m(3) - Icges(3,1) + Icges(3,2) + t56) * t13) * t13 + ((t32 * t55) + t18 + t19 + t40 + ((2 * pkin(2) + pkin(3)) * pkin(3)) + 0.2e1 * t30 * t6) * m(7); (-rSges(3,1) * t45 + Icges(3,5)) * t10 - (rSges(3,2) * t45 - Icges(3,6)) * t13 + (-(rSges(4,3) * m(4)) + t28 - t44) * t50; t56 + (t22 + t25) * m(3) + Icges(3,3); (-rSges(5,1) * t44 + Icges(5,5)) * t9 + (-rSges(5,2) * t44 + Icges(5,6)) * t12 + t28 * t54; ((t3 * t12 - t34) * t13 + (t12 * t47 + t3 * t9) * t10) * pkin(4); t57 + (t21 + t24) * m(5) + Icges(5,3); t30 * (t13 * pkin(4) + t43 + t6) * m(7) + t39; t50 * t53; t53 * t54; t39;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
