% Calculate joint inertia matrix for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1DE1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_inertiaJ_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1DE1_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1DE1_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:10
% EndTime: 2020-04-24 19:57:10
% DurationCPUTime: 0.21s
% Computational Cost: add. (308->39), mult. (425->63), div. (18->7), fcn. (110->4), ass. (0->32)
t46 = pkin(4) ^ 2;
t26 = cos(qJ(1));
t42 = pkin(2) * t26;
t40 = (-0.2e1 * t42 + pkin(1)) * pkin(1);
t43 = -pkin(3) + pkin(4);
t44 = -pkin(3) - pkin(4);
t45 = ((pkin(2) - t44) * (pkin(2) + t44) + t40) * ((pkin(2) - t43) * (pkin(2) + t43) + t40);
t25 = sin(qJ(1));
t41 = t25 * pkin(2);
t39 = pkin(3) ^ 2 - t46;
t31 = pkin(2) ^ 2;
t23 = t31 + t40;
t19 = 0.1e1 / t23;
t18 = t23 - t39;
t33 = sqrt(-t45);
t37 = -pkin(1) + t42;
t3 = 0.1e1 - pkin(1) * (t18 * t41 - t33 * t37) / t33 * t19;
t38 = t19 * t3 / pkin(3);
t34 = t23 ^ 2;
t22 = rSges(2,1) * t26 - rSges(2,2) * t25;
t21 = -rSges(2,1) * t25 - rSges(2,2) * t26;
t17 = t23 + t39;
t16 = -rSges(3,1) * t37 + rSges(3,2) * t41;
t15 = rSges(3,1) * t41 + rSges(3,2) * t37;
t14 = rSges(4,1) * t41 + rSges(4,2) * t37;
t13 = -rSges(4,1) * t37 + rSges(4,2) * t41;
t7 = (pkin(1) * t26 - pkin(2)) * t33 - pkin(1) * t25 * t17;
t5 = -t13 * t18 + t14 * t33;
t4 = t13 * t33 + t14 * t18;
t2 = t42 + (t15 * t33 + t16 * t17) * t38 / 0.2e1;
t1 = -t41 - (-t15 * t17 + t16 * t33) * t38 / 0.2e1;
t6 = [t3 ^ 2 * Icges(3,3) + m(2) * (t21 ^ 2 + t22 ^ 2) + Icges(2,3) + m(3) * (t1 ^ 2 + t2 ^ 2) + (-0.1e1 / t34 * Icges(4,3) + m(4) * (-t4 ^ 2 / 0.2e1 - t5 ^ 2 / 0.2e1) / t46 / t34 ^ 2 / 0.2e1) * t31 * t7 ^ 2 / t45;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t6(1);];
Mq = res;
