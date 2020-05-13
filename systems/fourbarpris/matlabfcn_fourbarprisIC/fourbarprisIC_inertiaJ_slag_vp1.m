% Calculate joint inertia matrix for
% fourbarprisIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisIC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_inertiaJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisIC_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisIC_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisIC_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:32
% EndTime: 2020-05-07 09:59:32
% DurationCPUTime: 0.05s
% Computational Cost: add. (40->23), mult. (81->37), div. (8->4), fcn. (54->4), ass. (0->15)
t23 = rSges(3,3) + qJ(2);
t21 = cos(qJ(1));
t17 = t21 ^ 2;
t18 = sin(qJ(3));
t19 = sin(qJ(1));
t20 = cos(qJ(3));
t22 = (-t20 * t18 - t21 * t19) / (pkin(3) + qJ(2)) / (-t20 ^ 2 + t17);
t14 = -t21 * rSges(2,1) + t19 * rSges(2,2);
t13 = -t20 * rSges(4,1) + t18 * rSges(4,2);
t12 = t19 * rSges(2,1) + t21 * rSges(2,2);
t11 = t18 * rSges(4,1) + t20 * rSges(4,2);
t9 = t18 * t21 - t19 * t20;
t8 = -t19 * rSges(3,1) - t23 * t21;
t7 = -t21 * rSges(3,1) + t23 * t19;
t1 = [m(3) * (t19 ^ 2 + t17) + 0.1e1 / pkin(2) ^ 2 / t9 ^ 2 * (m(4) * (t11 ^ 2 + t13 ^ 2) + Icges(4,3)) + (0.2e1 * m(3) * (-t19 * t8 - t21 * t7) + (Icges(2,3) + Icges(3,2) + m(2) * (t12 ^ 2 + t14 ^ 2) + m(3) * (t7 ^ 2 + t8 ^ 2)) * t22) * t22;];
Mq = t1;
