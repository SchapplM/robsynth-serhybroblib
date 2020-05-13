% Calculate joint inertia matrix for
% fourbarprisOL
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisOL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_inertiaJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisOL_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisOL_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:04
% EndTime: 2020-05-07 09:52:04
% DurationCPUTime: 0.07s
% Computational Cost: add. (21->16), mult. (50->27), div. (0->0), fcn. (24->4), ass. (0->12)
t11 = rSges(3,3) + qJ(2);
t10 = cos(qJ(1));
t9 = cos(qJ(3));
t8 = sin(qJ(1));
t7 = sin(qJ(3));
t6 = -t10 * rSges(2,1) + t8 * rSges(2,2);
t5 = -t9 * rSges(4,1) + t7 * rSges(4,2);
t4 = t8 * rSges(2,1) + t10 * rSges(2,2);
t3 = t7 * rSges(4,1) + t9 * rSges(4,2);
t2 = -t8 * rSges(3,1) - t11 * t10;
t1 = -t10 * rSges(3,1) + t11 * t8;
t12 = [Icges(2,3) + Icges(3,2) + m(2) * (t4 ^ 2 + t6 ^ 2) + m(3) * (t1 ^ 2 + t2 ^ 2); m(3) * (-t1 * t10 - t8 * t2); m(3) * (t10 ^ 2 + t8 ^ 2); 0; 0; m(4) * (t3 ^ 2 + t5 ^ 2) + Icges(4,3); 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t12(1), t12(2), t12(4), t12(7); t12(2), t12(3), t12(5), t12(8); t12(4), t12(5), t12(6), t12(9); t12(7), t12(8), t12(9), t12(10);];
Mq = res;
