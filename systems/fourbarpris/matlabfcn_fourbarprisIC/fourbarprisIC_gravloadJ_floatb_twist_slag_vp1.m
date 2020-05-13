% Calculate Gravitation load on the joints for
% fourbarprisIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbarprisIC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:32
% EndTime: 2020-05-07 09:59:32
% DurationCPUTime: 0.07s
% Computational Cost: add. (19->18), mult. (33->31), div. (4->4), fcn. (26->4), ass. (0->6)
t5 = rSges(3,3) + qJ(2);
t4 = cos(qJ(1));
t3 = cos(qJ(3));
t2 = sin(qJ(1));
t1 = sin(qJ(3));
t6 = [(-t3 * t1 - t4 * t2) / (pkin(3) + qJ(2)) / (-t3 ^ 2 + t4 ^ 2) * (-m(2) * (g(1) * (t2 * rSges(2,1) + t4 * rSges(2,2)) + g(2) * (-t4 * rSges(2,1) + t2 * rSges(2,2))) - m(3) * (g(1) * (-t4 * rSges(3,1) + t5 * t2) + g(2) * (-t2 * rSges(3,1) - t5 * t4))) - m(3) * (-g(1) * t4 - g(2) * t2) + 0.1e1 / pkin(2) / (t1 * t4 - t2 * t3) * m(4) * (g(1) * (t1 * rSges(4,1) + t3 * rSges(4,2)) + g(2) * (-t3 * rSges(4,1) + t1 * rSges(4,2)));];
taug = t6(:);
