% Implicit kinematic constraints of
% fourbar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = fourbar1IC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_kinconstr_impl_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_kinconstr_impl_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:48
% EndTime: 2020-04-24 20:15:48
% DurationCPUTime: 0.04s
% Computational Cost: add. (11->9), mult. (6->6), div. (0->0), fcn. (6->6), ass. (0->2)
t1 = qJ(1) + qJ(2);
t2 = [cos(qJ(3)) * pkin(4) + pkin(1) + pkin(3) * cos(t1) - cos(qJ(1)) * pkin(2); sin(qJ(3)) * pkin(4) + pkin(3) * sin(t1) - sin(qJ(1)) * pkin(2); pi + qJ(3) - qJ(4) - t1;];
h = t2;
