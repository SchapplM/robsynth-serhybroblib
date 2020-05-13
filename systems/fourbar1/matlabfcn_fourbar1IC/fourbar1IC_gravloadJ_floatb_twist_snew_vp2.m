% Calculate Gravitation load with newton euler on the joints for
% fourbar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:50
% EndTime: 2020-04-24 20:15:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->21), mult. (41->32), div. (4->3), fcn. (33->8), ass. (0->13)
t34 = -qJ(3) + qJ(1);
t33 = cos(qJ(1));
t32 = cos(qJ(2));
t31 = cos(qJ(3));
t30 = sin(qJ(1));
t29 = sin(qJ(2));
t28 = sin(qJ(3));
t27 = sin(qJ(2) + t34);
t25 = -t33 * g(1) - t30 * g(2);
t24 = t30 * g(1) - t33 * g(2);
t23 = t29 * mrSges(3,1) + t32 * mrSges(3,2) - mrSges(2,2);
t22 = m(3) * pkin(2) - t32 * mrSges(3,1) + t29 * mrSges(3,2) + mrSges(2,1);
t1 = [-(-t30 * t22 + t23 * t33) * g(1) - (t22 * t33 + t30 * t23) * g(2) + ((-pkin(3) * t27 + pkin(2) * sin(t34)) / pkin(3) * (mrSges(3,1) * (-t32 * t24 + t29 * t25) - mrSges(3,2) * (-t29 * t24 - t32 * t25)) + t29 * pkin(2) / pkin(4) * (-(-t28 * mrSges(4,1) - mrSges(4,2) * t31) * g(1) - (mrSges(4,1) * t31 - t28 * mrSges(4,2)) * g(2))) / t27;];
taug = t1(:);
