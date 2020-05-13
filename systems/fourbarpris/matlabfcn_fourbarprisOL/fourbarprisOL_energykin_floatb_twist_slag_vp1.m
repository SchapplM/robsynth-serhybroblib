% Calculate kinetic energy for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisOL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:02
% EndTime: 2020-05-07 09:52:03
% DurationCPUTime: 0.56s
% Computational Cost: add. (161->95), mult. (238->127), div. (0->0), fcn. (126->4), ass. (0->44)
t77 = Icges(2,4) - Icges(3,5);
t76 = Icges(2,1) + Icges(3,3);
t75 = Icges(3,1) + Icges(2,2);
t52 = sin(qJ(1));
t74 = t77 * t52;
t54 = cos(qJ(1));
t73 = t77 * t54;
t72 = t76 * t54 - t74;
t71 = t76 * t52 + t73;
t70 = -t75 * t54 - t74;
t69 = t75 * t52 - t73;
t68 = rSges(3,3) + qJ(2);
t67 = Icges(3,4) + Icges(2,6);
t66 = -Icges(2,5) + Icges(3,6);
t51 = sin(qJ(3));
t63 = Icges(4,4) * t51;
t53 = cos(qJ(3));
t62 = Icges(4,4) * t53;
t59 = V_base(6) * pkin(1) + V_base(2);
t56 = -t52 * rSges(3,1) - t68 * t54;
t55 = -rSges(3,1) * t54 + t68 * t52;
t49 = V_base(6) + qJD(1);
t48 = V_base(6) + qJD(3);
t47 = -rSges(2,1) * t54 + t52 * rSges(2,2);
t45 = -rSges(4,1) * t53 + rSges(4,2) * t51;
t44 = -t52 * rSges(2,1) - rSges(2,2) * t54;
t42 = -rSges(4,1) * t51 - rSges(4,2) * t53;
t37 = -Icges(4,1) * t53 + t63;
t36 = -Icges(4,1) * t51 - t62;
t31 = Icges(4,2) * t51 - t62;
t30 = -Icges(4,2) * t53 - t63;
t23 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t22 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t21 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t20 = V_base(5) * rSges(2,3) - t44 * t49 + V_base(1);
t19 = V_base(5) * rSges(4,3) - t42 * t48 + V_base(1);
t18 = -V_base(4) * rSges(4,3) + t45 * t48 + V_base(2);
t17 = -V_base(4) * rSges(2,3) + t47 * t49 + t59;
t16 = t42 * V_base(4) - t45 * V_base(5) + V_base(3);
t15 = t44 * V_base(4) + V_base(3) + (-pkin(1) - t47) * V_base(5);
t14 = -V_base(5) * rSges(3,2) - qJD(2) * t54 + t55 * t49 + V_base(1);
t13 = V_base(4) * rSges(3,2) - qJD(2) * t52 + t56 * t49 + t59;
t12 = V_base(3) - t55 * V_base(4) + (-pkin(1) - t56) * V_base(5);
t1 = m(1) * (t21 ^ 2 + t22 ^ 2 + t23 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t15 ^ 2 + t17 ^ 2 + t20 ^ 2) / 0.2e1 + m(3) * (t12 ^ 2 + t13 ^ 2 + t14 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((-Icges(4,5) * t51 - Icges(4,6) * t53) * V_base(5) + (-Icges(4,5) * t53 + Icges(4,6) * t51) * V_base(4) + Icges(4,3) * t48 / 0.2e1) * t48 + ((Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * t49 + (t66 * V_base(4) - t67 * V_base(5)) * t54 + (t66 * V_base(5) + t67 * V_base(4)) * t52) * t49 + ((t30 * t51 - t36 * t53 + t70 * t52 + t71 * t54) * V_base(5) + (t51 * t31 - t53 * t37 + t69 * t52 + t72 * t54) * V_base(4)) * V_base(4) / 0.2e1 + ((-t53 * t30 - t51 * t36 + t71 * t52 - t70 * t54) * V_base(5) + (-t31 * t53 - t37 * t51 + t72 * t52 - t69 * t54) * V_base(4)) * V_base(5) / 0.2e1;
T = t1;
