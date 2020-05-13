% Calculate joint inertia matrix for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fivebar1OL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1OL_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1OL_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:00
% EndTime: 2020-04-27 06:13:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (28->22), mult. (54->28), div. (0->0), fcn. (8->4), ass. (0->9)
t18 = rSges(5,1) ^ 2 + rSges(5,2) ^ 2;
t17 = rSges(3,1) * cos(qJ(2));
t16 = rSges(5,1) * cos(qJ(4));
t15 = rSges(3,2) * sin(qJ(2));
t14 = rSges(5,2) * sin(qJ(4));
t13 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t12 = t18 * m(5) + Icges(5,3);
t11 = t13 * m(3) + Icges(3,3);
t1 = [(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3) + ((pkin(2) + 0.2e1 * t15 - 0.2e1 * t17) * pkin(2) + t13) * m(3); (t15 - t17) * m(3) * pkin(2) + t11; t11; 0; 0; (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3) + Icges(5,3) + ((pkin(3) - 0.2e1 * t14 + 0.2e1 * t16) * pkin(3) + t18) * m(5); 0; 0; (-t14 + t16) * m(5) * pkin(3) + t12; t12; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
