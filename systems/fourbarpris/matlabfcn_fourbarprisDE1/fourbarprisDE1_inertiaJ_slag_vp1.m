% Calculate joint inertia matrix for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisDE1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_inertiaJ_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_inertiaJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE1_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE1_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE1_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:09:51
% EndTime: 2020-05-07 09:09:51
% DurationCPUTime: 0.26s
% Computational Cost: add. (142->52), mult. (205->71), div. (5->2), fcn. (6->2), ass. (0->27)
t11 = (qJ(1) + pkin(3));
t15 = rSges(3,3) ^ 2;
t19 = rSges(3,1) ^ 2;
t41 = (-t15 - t19);
t36 = -2 * rSges(3,3) * qJ(1) + t41;
t48 = -t15 / 0.2e1 - t19 / 0.2e1;
t47 = ((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4));
t45 = 2 * t11;
t37 = -pkin(2) + t11;
t38 = -pkin(2) - t11;
t35 = (pkin(1) + t38) * (pkin(1) + t37) * (pkin(1) - t37) * (pkin(1) - t38);
t43 = sqrt(-t35) * rSges(3,1);
t14 = qJ(1) ^ 2;
t42 = rSges(3,3) * t14;
t23 = pkin(1) ^ 2;
t39 = 2 * t23;
t34 = -(2 * t42) + t43;
t12 = rSges(2,2) ^ 2 * m(2);
t13 = m(2) * rSges(2,1) ^ 2;
t33 = (-t12 - t13 - Icges(3,2) - Icges(2,3));
t22 = (pkin(3) ^ 2);
t32 = 2 * pkin(3) * qJ(1) + t22 + t36;
t31 = t11 * t33;
t29 = t11 ^ 2;
t25 = pkin(2) ^ 2;
t21 = pkin(3) * t22;
t1 = [((-(t33 * t39) + (t31 * t45) + (-(t32 * t39) + (-t21 + ((-4 * t14 + t36) * pkin(3)) + t34 + ((-3 * t22 - 2 * t14 + t41) * qJ(1))) * t45) * m(3)) * t25 - 0.2e1 * ((t21 + ((2 * t14 + t36) * pkin(3)) + ((3 * t22 + t41) * qJ(1)) + t34) * m(3) + t31) * t11 * t23 + ((-t12 / 0.2e1 - t13 / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(2,3) / 0.2e1 - (2 * Icges(4,3)) - (2 * t47)) * pkin(3) - (qJ(1) * (4 * Icges(4,3) + 4 * t47 - t33)) / 0.2e1 + (t21 / 0.2e1 + (t14 + t48) * pkin(3) - t42 + t43 + (0.3e1 / 0.2e1 * t22 - (rSges(3,3) * pkin(3)) + t48) * qJ(1)) * m(3)) * t29 * t45 + (t23 ^ 2 + t25 ^ 2) * (t32 * m(3) + t33)) / t29 / t35;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t1(1);];
Mq = res;
